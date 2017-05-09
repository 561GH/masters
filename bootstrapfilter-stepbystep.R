####################################################################
# 
# BOOTSTRAP FILTER
#
#   y(t) = x(t)/(1+x(t)^2) + v(t)
#   x(t) = x(t-1) + w(t)
#
#   v(t) ~ N(0,1)
#   w(t) ~ N(0,0.5)
#   x(0) ~ N(1,10)
#
# Reference: Gordon, Salmond and Smith (1993)
#
####################################################################
#
# HEDIBERT FREITAS LOPES
# Associate Professor of Econometrics and Statistics
# The University of Chicago Booth School of Business
# 5807 South Woodlawn Avenue
# Chicago, Illinois, 60637
# Email : hlopes@ChicagoGSB.edu
# URL: http://faculty.chicagobooth.edu/hedibert.lopes/research/
#
####################################################################
q025 = function(x){quantile(x,0.025)}
q975 = function(x){quantile(x,0.975)}

####################################################################
# Step by step
####################################################################
set.seed(12254)
par(mfrow=c(2,5))

N = 1000
# Posterior at time t=0 p(x[0]|y[0])=N(1,10)
x = rnorm(N,1,sqrt(10))
hist(x,main="Pr(x[0]|y[0])")

# Obtain draws from prior p(x[1]|y[0])
x1 = x + rnorm(N,0,sqrt(0.5))
hist(x1,main="Pr(x[1]|y[0])")

# Let us take a look at the likelihood function
# p(y[1]|x[1])
y1  = 5
ths = seq(-30,30,length=1000)
plot(ths,dnorm(y1,ths/(1+ths^2),1),type="l",xlab="",ylab="")
title(paste("p(y[1]=",y1,"|x[1])",sep=""))

# Computing resampling weights
w = dnorm(y1,x1/(1+x1^2),1)

# Resample to obtain draws from p(x[1]|y[1])
k     = sample(1:N,size=N,replace=TRUE,prob=w)
x = x1[k]
hist(x,main="Pr(x[1]|y[1])")

# Obtain draws from prior p(x[2]|y[1])
x2 = x + rnorm(N,0,sqrt(0.5))
hist(x2,main="Pr(x[2]|y[1])")

# Let us take a look at the likelihood function
# p(y[2]|x[2])
y2  = -2
ths = seq(-30,30,length=1000)
plot(ths,dnorm(y2,ths/(1+ths^2),1),type="l",xlab="",ylab="")
title(paste("p(y[2]=",y2,"|x[2])",sep=""))

# Computing resampling weights
w = dnorm(y2,x2/(1+x2^2),1)

# Resample to obtain draws from p(x[2]|y[2])
k     = sample(1:N,size=N,replace=TRUE,prob=w)
x = x2[k]
hist(x,main="Pr(x[2]|y[2])")


# Obtain draws from prior p(x[3]|y[2])
x3 = x + rnorm(N,0,sqrt(0.5))
hist(x3,main="Pr(x[3]|y[2])")

# Let us take a look at the likelihood function
# p(y[3]|x[3])
y3  = 20
ths = seq(-30,30,length=1000)
plot(ths,dnorm(y3,ths/(1+ths^2),1),type="l",xlab="",ylab="")
title(paste("p(y[3]=",y3,"|x[3])",sep=""))

# Computing resampling weights
w = dnorm(y3,x3/(1+x3^2),1)

# Resample to obtain draws from p(x[3]|y[3])
k     = sample(1:N,size=N,replace=TRUE,prob=w)
x = x3[k]
hist(x,main="Pr(x[3]|y[3])")

######################################################
# Studying particle size and associated MC error
######################################################
set.seed(1235)
n = 100
x = rep(0,n)
y     = rep(0,n)
x[1] = 2
y[1] = rnorm(1,x[1]/(1+x[1]^2),1)
for (t in 2:n){
  x[t] = x[t-1] + rnorm(1,0,sqrt(0.5))
  y[t]     = rnorm(1,x[t]/(1+x[t]^2),1)
}
par(mfrow=c(1,2))
ts.plot(x)
ts.plot(y)
xtrue = x

set.seed(923579)
nrep =  10
L    = -10
U    =  20
Ns   = c(1000,2000,5000,10000)
par(mfrow=c(2,2))
for (N in Ns){
  for (rep in 1:nrep){
    xs = matrix(0,n,N)
    x = rnorm(N,1,sqrt(10))
    for (t in 1:n){
      x1 = x + rnorm(N,0,sqrt(0.5))
      w      = dnorm(y[t],x1/(1+x1^2),1)
      k      = sample(1:N,size=N,replace=TRUE,prob=w)
      x = x1[k]
      xs[t,] = x
    }
    mx = apply(xs,1,median)
    lx = apply(xs,1,q025)
    ux = apply(xs,1,q975)
    if (rep==1){
      plot(mx,type="l",ylim=c(L,U),xlab="Time",col=rep)
      title(paste("N=",N,sep=""))
      lines(lx,col=rep)
      lines(ux,col=rep)
    }else{
      lines(lx,col=rep)
      lines(mx,col=rep)
      lines(ux,col=rep)
    }
  }
}


