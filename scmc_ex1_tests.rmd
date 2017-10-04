# More 1-dimensional tests

```{r log(20 * x + 1)}
ytrue <- function(x) {log(20 * x + 1)}  # x > -1/20; monotone increasing function

# x = inputs where fcn evaluated
# xstar = prediction set inputs
# xprime = derivative set inputs
# true/known model values
given <- list(x = cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)),
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)) -
                mean(ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))) )

initial.values <- list(l = 0.5, sig2 = 10,
                       ystar = ytrue(given$xstar) - mean(ytrue(given$x)),
                       yprime = 20 / (20 * given$xprime) )
```

```{r 20*log(20*x+1)}
ytrue <- function(x) {20*log(20*x+1)}

# x = inputs where fcn evaluated
# xstar = prediction set inputs
# xprime = derivative set inputs
# true/known model values
given <- list(x = cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)),
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)) -
                mean(ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))) )

plot(given$x, given$y)

initial.values <- list(l = 0.5, sig2 = 10,
                       ystar = ytrue(given$xstar) - mean(ytrue(given$x)),
                       yprime = 20^2 / (20 * given$xprime) )

# seems like all particles die by the end so it's just selecting the same single particle eventually
# i.e. at some point algorithm dies since there's really only 1 particle
# this doesn't seem to be the problem if I follow the paper
```

```{r sin(x^10)}
ytrue <- function(x) {sin(x^10)}
yprimetrue <- function(x) {cos( x^10 ) * 10 * x^9 }

x <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.9, 1)
given <- list(x = cbind(x),
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(x) -
                mean(ytrue(x)) )
rm(x)

plot(given$x, given$y)


initial.values <- list(l = 0.5, sig2 = 10,
                       ystar = ytrue(given$xstar) - mean(ytrue(given$x)),
                       yprime = yprimetrue(given$xprime) )

# happens here too
```

```{r }
ytrue <- function(x) {x + 0.2*cos(x) + sin(x)}
yprimetrue <- function(x) {1 - 0.2*sin(x) + cos(x) }

x <- c(seq(0, 15, length.out = 10), seq(20, 25, length.out = 5))
given <- list(x = cbind(x) / 25,
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(x) -
                mean(ytrue(x)) )
rm(x)

plot(given$x, given$y, type = "l")
points(given$x, given$y)


initial.values <- list(l = 0.5, sig2 = 10,
                       ystar = ytrue(given$xstar) - mean(ytrue(given$x)),
                       yprime = yprimetrue(given$xprime) )

# happens here too
```


