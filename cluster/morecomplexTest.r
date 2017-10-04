# Try installing packages, etc...

getwd() # Just to see :)

install.packages('binom',lib='Packages',repos = "http://cran.us.r-project.org")
library(binom,lib.loc="Packages")

# Do some basic work
binom.confint(100,200)
