\name{spati2}
\alias{spati2}
\title{Heights and diameters of Scots pine trees in Ilomantsi, Finland.}
\description{A dataset of Scots pine tree heights and diameters. 
The trees were collected on 56 fixed-area sample plots. This is a subset of the larger data set \code{spati}.}
\usage{data(spati2)}
\format{
  A data frame with 1678 observations on the following 3 variables.
  \describe{
    \item{\code{plot}}{A unique sample plot id.}
    \item{\code{d}}{Tree diameter at breast height (1.3 meters above the ground) in cm.}
    \item{\code{h}}{Tree height, m.}
    \item{\code{n}}{The total number of trees on the plot.}
    \item{\code{dvar}}{The variance of tree diameters on the plot.}
    \item{\code{dmean}}{The mean of tree diameters on the plot.}
  }
}

\author{
The data were collected by Timo Pukkala.
}

\references{
Pukkala, T. 1989. Prediction of tree diameter and height in a Scots pine stand as a
function of the spatial pattern of trees. Silva Fennica 23(2): 83-99.
}

\examples{
data(spati2)
fithd(spati2$d,spati2$h,spati2$plot)
}

\keyword{ datasets }
