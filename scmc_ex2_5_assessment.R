library(sm)

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
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(20)

sm.density(particles.l[,1,1], col = cols[1])
abline(v = mean(particles.l[,1,1]), col = cols[1])

for (t in 2:M) {
  sm.density(particles.l[,t,1], col = cols[t], add = TRUE)
  abline(v = mean(particles.l[,t,1]), col = cols[t])
}

sm.density(particles.l[,1,1], col = cols[1], add = TRUE)
abline(v = mean(particles.l[,1,1]), col = cols[1])

## l2
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(20)

sm.density(particles.l[,1,2], col = cols[1])
abline(v = mean(particles.l[,1,2]), col = cols[1])

for (t in 2:M) {
  sm.density(particles.l[,t,2], col = cols[t], add = TRUE)
  abline(v = mean(particles.l[,t,2]), col = cols[t])
}

sm.density(particles.l[,1,2], col = cols[1], add = TRUE)
abline(v = mean(particles.l[,1,2]), col = cols[1])


# assess sig2 #####################
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(20)

sm.density(particles.sig2[,1], col = cols[1])
abline(v = mean(particles.sig2[,1]), col = cols[1])

for (t in 2:M) {
  sm.density(particles.sig2[,t], col = cols[t], add = TRUE)
  abline(v = mean(particles.sig2[,t]), col = cols[t])
}

sm.density(particles.sig2[,1], col = cols[1], add = TRUE)
abline(v = mean(particles.sig2[,1]), col = cols[1])

# acceptance rates? #################
all.rates <- NULL
for (t in 1:M) {
  all.rates <- rbind(all.rates,
                     (acceptances[[t]])[N,])
}
all.rates.frame <- cbind(time = 2:M, all.rates / N)

# FIX
