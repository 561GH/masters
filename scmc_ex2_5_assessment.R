# assess l #####################
# sm.density(particles.l[,M,])
# sm.density(particles.l[,M,], display = "slice")
# sm.density(particles.l[,M,], display = "image")

for (t in 1:M) {
  pdf(file = paste("time_", t, ".pdf", sep = ""))
  sm.density(particles.l[,t,], display = "image")
  dev.off()
}

## l1
sm.density(particles.l[,1,1], col = "red")
abline(v = mean(particles.l[,1,1]), col = "red")

sm.density(particles.l[,1,1], col = "red")
for (t in 2:(M-1)) {
  sm.density(particles.l[,t,1], add = TRUE)
  abline(v = mean(particles.l[,t,1]))
}

sm.density(particles.l[,1,1], col = "red")
abline(v = mean(particles.l[,1,1]), col = "red")
sm.density(particles.l[,M,1], col = "green", add = TRUE)
abline(v = mean(particles.l[,M,1]), col = "green")

## l2
sm.density(particles.l[,1,2], col = "red")
abline(v = mean(particles.l[,1,2]), col = "red")

sm.density(particles.l[,1,2], col = "red")
for (t in 2:(M-1)) {
  sm.density(particles.l[,t,2], add = TRUE)
  abline(v = mean(particles.l[,t,2]))
}

sm.density(particles.l[,1,2], col = "red")
abline(v = mean(particles.l[,1,2]), col = "red")
sm.density(particles.l[,M,2], col = "green", add = TRUE)
abline(v = mean(particles.l[,M,2]), col = "green")


# assess sig2 #####################
sm.density(particles.sig2[,1], col = "red")

for (t in 2:(M-1)) {
  sm.density(particles.sig2[,t], add = TRUE)
  abline(v = mean(particles.sig2[,t]))
}

sm.density(particles.sig2[,1], col = "red", add = TRUE)
abline(v = mean(particles.sig2[,1]), col = "red")

sm.density(particles.sig2[,M], col = "green", add = TRUE)
abline(v = mean(particles.sig2[,M]), col = "green")
