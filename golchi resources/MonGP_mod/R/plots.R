
par_density = function(obj, par) {
  if (par == 'l') {
    df0 = as.data.frame(obj$l)
    df0 = melt(df0)
    names(df0) = c('step', 'l')
    L = length(levels(df0$step))
    levels(df0$step) = 1:L
    select_step = floor(seq(1, L, length = 5))
    df0 = df0 %>%
      filter(step %in% select_step)
    p = ggplot(df0, aes(x = l, fill = step)) +
      geom_density(alpha = .5, color = 'grey') + scale_fill_brewer()
  }
  if (par == 'sig2') {
    df0 = as.data.frame(obj$sig2)
    df0 = melt(df0)
    names(df0) = c('step', 'sig2')
    L = length(levels(df0$step))
    levels(df0$step) = 1:L
    select_step = floor(seq(1, L, length = 5))
    df0 = df0 %>%
      filter(step %in% select_step)
    p = ggplot(df0, aes(x = sig2, fill = step)) +
      geom_density(alpha = .5, color = 'grey') + scale_fill_brewer()
  }
  return(p)
}

#' Plot sample paths
#'
#'
#' @param obj output of SMC_MonGP
#' @export
#' @return ggplott object
#'

sample_path = function(obj) {
  L = dim(obj$l)[2]
  N = dim(obj$l)[1]
  dimnames(obj$y) = list(sample = 1:N, step = 1:L, point = obj$newx)
  df0 <- as.data.frame.table(obj$y)
  names(df0)[4] = 'y'
  select_step = floor(seq(1, L, length = 5))
  df0 = df0 %>%
    filter(step %in% select_step) %>%
    group_by(point, step) %>%
    summarise(mean = mean(y), lower = quantile(y, .025), upper = quantile(y, .975))
  p = ggplot(df0, aes(x = as.numeric(as.character(point)), y = mean, fill = step)) + geom_line() +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) + scale_fill_brewer() + xlab('x') + ylab('y')
  return(p)
}


