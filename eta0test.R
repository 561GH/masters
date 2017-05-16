post.args <- list(l = 5,
                  sig2 = 1,
                  ystar = ytrue(xstar),
                  yprime = 20 / (20 * xprime),
                  y = y,
                  xstar = cbind(xstar), xprime = cbind(xprime), x = cbind(x) )

do.call(posterior, post.args)
