\name{plants2}
\alias{plants2}
\docType{data}
\title{
Sapling counts from sample plots of sapling stands in Finland.
}
\description{Independent Norway Spruce regeneration establishment data.}
\usage{data(plants2)}
\format{
  A data frame with 123 observations on the following 8 variables.
  \describe{
    \item{\code{planted}}{The number of planted spruce saplings on the plot}
    \item{\code{pines}}{The number of natural pine saplings}
    \item{\code{spruces}}{The number of natural spruce saplings}
    \item{\code{birches}}{The number natural birch  saplings}
    \item{\code{othersp}}{The number of natural saplings of other species}
    \item{\code{hcrop}}{The mean height of crop species}
    \item{\code{hdecid}}{The mean height of deciduous tree species}
    \item{\code{sitetype}}{Site fertility class (small number indicates more fertile site)}
  }
}
\details{
The data are collected from 123 fixed-area sample plots with similar age of planted spruce saplings. The number of saplings per species and the height of 
crop species (spruce and pine) and competing vegetation (birch and other broadleaved trees) has been recorded for all plots. 
The data includes one plot per forest stand.
}

\references{
Miina, J. and Saksa, T. 2006. Predicting regeneration establishment in Norway spruce plantations using a multivariate multilevel model. New Forests 32: 265-283.
}

\examples{
data(plants2)

glm1 <- glm(spruces ~ hdecid, family=quasipoisson(), data=plants2)
}
\keyword{datasets}
