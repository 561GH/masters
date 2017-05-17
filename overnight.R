# Running code overnight

install.packages(c("foreach", "mvtnorm"))
install.packages(file.choose(), repos = NULL, type = "source", dependencies = TRUE)


saveRDS(chains, file = "C:\Users\gghsu\Desktop\eta0chains.rds")
saveRDS(out, file = "C:\Users\gghsu\Desktop\eta0out.rds")
