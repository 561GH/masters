# SCMC toy example
## response y must be scalar

# run: scmc_ex1_initialise
# + build package
set.seed(761)

M <- 20  # total time
# tauM <- 10e-6
# taus <- seq(0, tauM, length.out = M+1)  # TODO: sequence of taus???
# taus <- taus[2:length(taus)]

# shirin's tau sequence i.e. reciprocal of her nuseq
taus <- 1 / c(Inf, ( seq(2, .1, length.out = M-1) )^5)

N <- 100  # particles desired for SCMC
burn <- 300  # burn in for initialising particles

################################################################################
# 1. initialise (ith row = ith particle) #######################################
################################################################################
# dim1 = particle, dim2 = time, dim3 = input index
particles.l <- matrix(NA, nrow = N, ncol = M)
particles.sig2 <- matrix(NA, nrow = N, ncol = M)
particles.ystar <- array(NA, dim = c(N, M, nrow(given$xstar)))
particles.yprime <- array(NA, dim = c(N, M, nrow(given$xprime)))

init <- eta0(eta.init = list(l = 0.5, sig2 = 0.1,
                             ystar = ytrue(given$xstar) - mean(ytrue(given$x)) +
                               rnorm(nrow(given$xstar), sd = 0.1),  # so init is not so good
                             yprime = 20 / (20 * given$xprime) ),
             given = given,  # data, locations, obs, etc.)
             N = burn + N, # particles
             v1 = 0.05, # step size for l proposal
             v2 = 0.01) # step size for sig2 proposal

particles.l[,1] <- (init$chain.l)[-(1:burn)]
particles.sig2[,1] <- (init$chain.sig2)[-(1:burn)]
particles.ystar[,1,] <- (init$chain.ystar)[-(1:burn),]
particles.yprime[,1,] <- (init$chain.yprime)[-(1:burn),]

cat("Acceptance rates from initialisation: \n",
    "\t l: ", sum(init$accepted.l[-(1:burn)]) / N, "\n",
    "\t sig2: ", sum(init$accepted.sig2[-(1:burn)]) / N, "\n \n")

################################################################################
# 2. weights at time t = 1 #####################################################
################################################################################
## matrix of weights: jth column = weights for t-th iteration
## dim1 = particle, dim2 = t-th iteration weight
W <- matrix(NA, nrow = N, ncol = M)
W[,1] <- 1/N

################################################################################
# 3. looping through t = 1, ..., M-1 ##########################################
################################################################################
step.l <- 0.05  # step size for l
step.sig2 <- 0.01  # step size for sig2
step.ystaryprime <- 1  # tune this

# TODO: should let step size vary to provide decent acceptance rates (during burn in?)

# ==============================================================================
## despite notation in the algorithm, "eta0" will really count as the first
## iteration of SCMC just to make the code cleaner
## => hence time t goes from 2 to M

for (t in 2:M) {

  cat("\n Starting time step:", t, "out of", M, "\n")

  particles.l[,t] <- particles.l[,t-1]
  particles.sig2[,t] <- particles.sig2[,t-1]
  particles.ystar[,t,] <- particles.ystar[,t-1,]
  particles.yprime[,t,] <- particles.yprime[,t-1,]

  ## WEIGHT CALCULATION ========================================================
  ### dim1 = particle, dim2 = t-th iteration weight

  ### calculate weight for each particle

  # TODO: debugging the dying weights issue
  # for (i in 1:N) {
  #   x <- particles.yprime[i,t,]
  #   num <- sum( log( pnorm( x * taus[t] )))
  #   den <- sum( log( pnorm( x * taus[t-1] )))
  #   cat("num:", num, "; den:", den, "\n")
  # }

  weight.particle <- function(x) {
    # x = derivative values for all inputs for a particle
    # i.e. vector of y'(x_d) for all d
    # print(length(x))  # make sure this is the number of derivative inputs

    num <- sum( log( pnorm( x * taus[t] )))
    den <- sum( log( pnorm( x * taus[t-1] )))
    return( exp(num - den) )
    #return( (num - den) )
  }

  w.tildes <- apply(particles.yprime[,t,], MARGIN = 1, FUN = weight.particle)
  W[,t] <- W[,t-1] * w.tildes

  W[,t] <- W[,t] / sum(W[,t])  # normalise weights

  ## RESAMPLING ================================================================
  ### effective sample size
  ESSt <- 1 / sum( (W[,t])^2 )

  if (ESSt <= N/2) {

    resample <- sample(1:N, replace = TRUE, prob = W[,t])

    particles.l[,t] <- particles.l[resample,t]
    particles.sig2[,t] <- particles.sig2[resample,t]
    particles.ystar[,t,] <- particles.ystar[resample,t,]
    particles.yprime[,t,] <- particles.yprime[resample,t,]

    W[,t] <- 1/N

  }

  ## SAMPLING ==================================================================
  ## the t-th time step particles (not t+1 as the alg says)

  ### Track acceptance rate ====================================================
  accepted.l <- accepted.sig2 <- accepted.ystaryprime <- 0

  ### Sample and update chains =================================================
  # TODO: fill this in after update.all done (acceptance rates ??? ignore for now)
  # TODO: adaptively chance step sizes

  particles.new <- update.all(particles.l.old = particles.l[,t-1],
                              particles.sig2.old = particles.sig2[,t-1],
                              particles.ystar.old = particles.ystar[,t-1,],
                              particles.yprime.old = particles.yprime[,t-1,])

  acceptances <- particles.new$acceptances  # N by 3 matrix of acceptances

  particles.l[,t] <- particles.new$l
  particles.sig2[,t] <- particles.new$sig2
  particles.ystar[,t,] <- particles.new$ystar
  particles.yprime[,t,] <- particles.new$yprime

}
