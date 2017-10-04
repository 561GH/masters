# Install these packages and load them
install.packages(
  c("ggplot2",
    "dplyr",
    "dtplyr",
    "data.table",
    "plotly",
    "plot3D",
    "boot",
    "gridExtra",
    "binom"),lib='Packages',repos = "http://cran.us.r-project.org")


#library(spatstat,lib.loc="Packages")
library(ggplot2,lib.loc="Packages")
library(dplyr,lib.loc="Packages")
library(dtplyr,lib.loc="Packages")
library(data.table,lib.loc="Packages")
library(plotly,lib.loc="Packages")
library(plot3D,lib.loc="Packages")
library(boot,lib.loc="Packages")
library(gridExtra,lib.loc="Packages")
library(binom,lib.loc="Packages")