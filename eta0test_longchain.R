# Why is eta0 sucking so much

longeta0 <- readRDS("/Users/ggh/Desktop/eta0chainslong.rds")

chain.ystar <- longeta0$chain.ystar
chain.yprime <- longeta0$chain.yprime
plot(longeta0$chain.l, type = "l")  # poor mixing for l??
plot(longeta0$chain.sig2[1000:15000], type = "l")

###
star <- function(i) {
  plot(chain.ystar[i,], type = "l")
  abline(h = chain.ystar[i,1], col = "red")
}

all.equal(chain.ystar[10,], chain.ystar[20,])
star(25)

###
prime <- function(i) {
  plot(chain.yprime[i,], type = "l")
  abline(h = chain.yprime[i,1], col = "red")
}

prime(1)

par(mfrow = c(3, 3))
for (i in 2:10) {prime(i)}
par(mfrow = c(1, 1))

###
tmp <- chain.ystar[,12000:15000]
submean <- function(v) {v - mean(v)}
tmp <- apply(tmp, MARGIN = 2, submean)

star2 <- function(i) {
  plot(tmp[i,], type = "l")
  abline(h = chain.ystar[i,1], col = "red")
}

star2(7)
