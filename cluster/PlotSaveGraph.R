# Try to plot a graph and save it

t <- seq(0,2,.0001)
A <- 3.1
a <- 56
d <- 7.2
B <- 3.7
b <- 29

x <- A*sin(a*t+d)
y <- B*sin(b*t)

# Do plot
jpeg('CoolPlot1.jpg')
plot(x=x,y=y,type='l')
dev.off()
