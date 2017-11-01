# assess GP parameter posteriors ##################
```{r}
library(sm)

# assess l #####################
# sm.density(particles.l[,M,])
# sm.density(particles.l[,M,], display = "slice")
# sm.density(particles.l[,M,], display = "image")

# for (t in 1:M) {
#   pdf(file = paste("time_", t, ".pdf", sep = ""))
#   sm.density(particles.l[,t,], display = "image")
#   dev.off()
# }

par(mfrow = c(2, 2))

## l1
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(M)

sm.density(particles.l[,1,1], col = cols[1], xlab = "l1")
abline(v = mean(particles.l[,1,1]), col = cols[1])

for (t in 2:M) {
  sm.density(particles.l[,t,1], col = cols[t], add = TRUE)
  abline(v = mean(particles.l[,t,1]), col = cols[t])
}

sm.density(particles.l[,1,1], col = cols[1], add = TRUE)
abline(v = mean(particles.l[,1,1]), col = cols[1])

## l2
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(M)

sm.density(particles.l[,1,2], col = cols[1], xlab = "l2")
abline(v = mean(particles.l[,1,2]), col = cols[1])

for (t in 2:M) {
  sm.density(particles.l[,t,2], col = cols[t], add = TRUE)
  abline(v = mean(particles.l[,t,2]), col = cols[t])
}

sm.density(particles.l[,1,2], col = cols[1], add = TRUE)
abline(v = mean(particles.l[,1,2]), col = cols[1])


# assess sig2 #####################
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(M)

#particles.sig2 <- sqrt(particles.sig2)

sm.density(particles.sig2[,1], col = cols[1], xlab = "sig2")
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
all.rates.frame <- data.frame(cbind(time = 2:M, all.rates / N))

plot(x = all.rates.frame$time, y = all.rates.frame$accepted.l,
     type = "l", col = "red",
     xlab = "time", ylab = "rate",
     ylim = c(min(all.rates.frame[,-1]), max(all.rates.frame[,-1])),
     main = paste("acceptance rates; time 2 to", M))

points(x = all.rates.frame$time, y = all.rates.frame$accepted.sig2,
       type = "l", col = "blue")

points(x = all.rates.frame$time, y = all.rates.frame$accepted.ystaryprime,
       type = "l", col = "darkgreen")

legend(x = 1.5, y = 0.25,
       legend = c("l", "sig2", "ystaryprime"),
       lty = "solid", lwd = 1, cex = 0.5,
       col = c("red", "blue", "darkgreen"),
       ncol = 3)

par(mfrow = c(1, 1))
```

# assess ystar posterior #################
```{r}
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
#colfunc <- colorRampPalette(c("white", "black"))
cols <- colfunc(M)

par(mfrow = c(3, 2))

for (i in 1:dim(particles.ystar)[3]) {

  parts.i <- particles.ystar[,,i]

  # for ylab
  ymax <- sm.density(parts.i[,M], col = cols[M], display = "none")

  sm.density(parts.i[,1], col = cols[1],
             ylim = c(0, max(ymax$estimate)),
             xlim = c(true.values$ystar[i] - 1.5, true.values$ystar[i] + 1.5),
             xlab = paste("ystar", i))
  abline(v = mean(parts.i[,1]), col = cols[1])

  for (t in 2:M) {
    sm.density(parts.i[,t], col = cols[t], add = TRUE)
    abline(v = mean(parts.i[,t]), col = cols[t])
  }

  abline(v = true.values$ystar[i], col = "black", lwd = 2)

  # sm.density(parts.i[,1], col = cols[1], add = TRUE)
  # abline(v = mean(parts.i[,1]), col = cols[1])

}

par(mfrow = c(1, 1))
```

# assess sig2 #################
```{r}
# tmp <- debug[[M]]
# head(tmp)
#
# prob.sig2 <- replace(tmp[,1], list = which(tmp[,1] > 1), values = 1)
# sig2new <- tmp[,2]
#

frame.prob.sig2 <- matrix(NA, nrow = N, ncol = M - 1)
frame.sig2new   <- matrix(NA, nrow = N, ncol = M - 1)

for (t in 2:M) {
  tmp <- debug[[t]]
  frame.prob.sig2[,t-1] <- replace(tmp[,1], list = which(tmp[,1] > 1), values = 1)
  frame.sig2new[,t-1] <- tmp[,2]
}

par(mfrow = c(1,1))
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(M-1)

plot(x = frame.sig2new[,1], y = frame.prob.sig2[,1],
     pch = 16, col = cols[1],
     xlim = c(min(frame.sig2new), max(frame.sig2new)),
     ylim = c(min(frame.prob.sig2), max(frame.prob.sig2)))

for (t in 2:ncol(frame.sig2new)) {
  points(x = frame.sig2new[,t], y = frame.prob.sig2[,t],
         pch = 16, col = cols[t])
}
```

# assess ystaryprime #################
```{r}
# tmp <- debug[[M]]
# head(tmp)
# 
# prob.ystarynew <- replace(tmp[,3], list = which(tmp[,3] > 1), values = 1)
# ystar2new <- tmp[,5]


frame.prob.ystarynew <- matrix(NA, nrow = N, ncol = M - 1)
frame.ystar2new   <- matrix(NA, nrow = N, ncol = M - 1)

for (t in 2:M) {
  tmp <- debug[[t]]
  frame.prob.ystarynew[,t-1] <- replace(tmp[,3], list = which(tmp[,3] > 1), values = 1)
  frame.ystar2new[,t-1] <- tmp[,5]
}

frame.prob.ystarynew <- replace(frame.prob.ystarynew, 
                                list = which(is.na(frame.prob.ystarynew)), 
                                values = 0)

par(mfrow = c(1,1))
colfunc <- colorRampPalette(c("red", "yellow", "springgreen", "royalblue"))
cols <- colfunc(M-1)

plot(x = frame.ystar2new[,1], y = frame.prob.ystarynew[,1],
     pch = 16, col = cols[1],
     xlim = c(min(frame.ystar2new), max(frame.ystar2new)),
     ylim = c(min(frame.prob.ystarynew), max(frame.prob.ystarynew)))

for (t in 2:ncol(frame.ystar2new)) {
  points(x = frame.ystar2new[,t], y = frame.prob.ystarynew[,t],
         pch = 16, col = cols[t])
}

abline(v = true.values$ystar[2], cex = 6)
```
