plotSamplePaths <- function(given,
                        particles.l,
                        particles.ystar) {
  L = dim(particles.l)[2]
  N = dim(particles.l)[1]
  dimnames(particles.ystar) = list(sample = 1:N, step = 1:L, point = given$xstar)
  df0 <- as.data.frame.table(particles.ystar)
  names(df0)[4] = "y"
  select_step = floor(seq(1, L, length = 5))
  df0 = df0 %>% filter(step %in% select_step) %>% group_by(point,
                                                           step) %>% summarise(mean = mean(y), lower = quantile(y,
                                                                                                                0.025), upper = quantile(y, 0.975))
  p = ggplot(df0, aes(x = as.numeric(as.character(point)),
                      y = mean, fill = step)) + geom_line() + geom_ribbon(aes(ymin = lower,
                                                                              ymax = upper), alpha = 0.2) + scale_fill_brewer() + xlab("x") +
    ylab("y")
  return(p)
}
