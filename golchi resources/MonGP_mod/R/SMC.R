#' SMC sampling to generate monotone GP sample paths
#'
#'
#' @param x Evaluation sites
#' @param y_obs Evaluations
#' @param xd virtual derivative sites
#' @param newx prediction sites
#' @param N Number of particles
#' @param L Number of densities in the sequence
#' @param initial_values a named list of initial values for the covariance parameters \code{l}
#' and \code{sig2}
#' @export
#' @return a named list: samples drawn from the sequence of posteriors for each parameter
#' @examples
#' x <- as.matrix(c(seq(0,.4,.1), seq(.9,1,.1)))
#' fnct <- function(x) {
#'   out <- log(20 * x + 1)
#'   return(out)
#' }
#' data <- fnct(x)
#' mu <- mean(data)
#' se <- as.numeric(sqrt(var(data)))
#' y_obs <- (data - mu)/se
#' newx <- as.matrix(seq(.01,.99,length=50))
#' xd <- as.matrix(seq(.42,.87,length=10))
#' initial_values = list(l = .5, sig2 = 10)
#' obj = SMC_MonGP(x, y_obs, xd, newx, N = 2000, L = 10, initial_values)
#' sample_path(obj)


SMC_MonGP = function(x, y_obs, xd, newx, N = 100, L = 10, initial_values){
  debug <- vector("list", L)
  weights <- matrix(NA, N, L)
  nuseq <- c(Inf, (seq(2, .1, length = L-1))^5)
  n = nrow(newx)
  m = nrow(xd)
  ESS <- c()
  lt <- array(dim = c(N, L))
  sig2t <- array(dim = c(N, L))
  yt <- array(dim = c(N, L, n))
  ydt <- array(dim = c(N, L, m))
  sample0 = MCMC_unconstrained(x, y_obs, xd, newx, burn = 200, N = 200 + N, initial_values)
  lt[,1] = sample0$l
  sig2t[,1] = sample0$sig2
  yt[,1,] = sample0$y
  ydt[,1,] = sample0$yd
  for (t in 2:L) {
    yt[,t,] <- yt[,t-1,]
    ydt[,t,] <- ydt[,t-1,]
    lt[,t] <- lt[,t-1]
    sig2t[,t] <- sig2t[,t-1]
    w = c()
    for (i in 1:N) {
      num <- sum(log(pnorm(ydt[i,t,]/nuseq[t])))
      den <- sum(log(pnorm(ydt[i,t,]/nuseq[t-1])))
      w[i] <- num - den
    }
    w = exp(w) / sum(exp(w))
    weights[,t] <- w
    ESS[t] <- 1/sum(w^2)
    index <- sample(1:N, N, replace = T, prob = w)
    yt[,t,] <- yt[index,t,]
    ydt[,t,] <- ydt[index,t,]
    lt[,t] <- lt[index,t]
    sig2t[,t] <- sig2t[index,t]
    #sampling new particles
    s0 = sampling(lt[,t], sig2t[,t], yt[,t,], ydt[,t,], M = 1, nu = nuseq[t])
    yt[,t,] <- s0$y
    ydt[,t,] <- s0$yd
    lt[,t] <- s0$l
    sig2t[,t] <- s0$sig2
    debug[[t]] <- s0$debug
  }
  out <- list(l = lt, sig2 = sig2t, y = yt, yd = ydt, x = x, newx = newx, y_obs = y_obs,
              debug = debug, weights = weights)
}

