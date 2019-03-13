\name{plants}
\alias{plants}
\docType{data}
\title{
Sapling counts from sample plots of sapling stands in Finland.
}
\description{Grouped Norway Spruce regeneration establishment data.}
\usage{data(plants)}
\format{
  A data frame with 1926 observations (fixed-area sample plots) from a total of 123 forest stands, with the following 8 variables.
  \describe{
    \item{\code{spruces}}{The number of spruce saplings (both planted and natural)}
    \item{\code{stand}}{The stand id)}
    \item{\code{hdecid}}{The mean height of deciduous tree species}
    \item{\code{prepar}}{Site preparation method, categorical with 4 levels}
    \item{\code{stones}}{Binary indicator for stoniness}
    \item{\code{wet}}{Binary indicator for wetness}
  }
}
\details{
The data are collected from 123 fixed-area sample plots with similar age of planted spruce saplings. The variables have been measured on fixed-area plots.
}

\references{
Miina, J. and Saksa, T. 2006. Predicting regeneration establishment in Norway spruce plantations using a multivariate multilevel model. New Forests 32: 265-283.
}

\examples{
data(plants)
library(lme4)
glmm1<-glmer(spruces ~ (1|stand)+hdecid+as.factor(prepar)+as.factor(stones)+as.factor(wet), 
             family=poisson(), 
			 data=plants)
}
\keyword{datasets}
