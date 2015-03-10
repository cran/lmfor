\name{spati}
\alias{spati}
\title{Growth of Scots pine trees in Ilomantsi, Finland}
\description{A dataset of Scots pine tree growth.}
\usage{spati}
\format{
\describe{
  \item{x,y}{the x- and y- coordinates of the midpoints}
  \item{r}{the radii of the circles}
  \item{col, lty, lwd}{the draving color, line type and line width of the perimeter line}
  \item{grayfill}{Logical- whether gray shading will be used for the circles or not}
  }
}
\value{
This function is used for its side effects.
}

\author{
Timo Pukkala 
}

\references{
Pukkala, T. 1989. Prediction of tree diameter and height in a Scots pine stand as a
function of the spatial pattern of trees. Silva Fennica 23(2): 83-99.
}

\examples{
plot(0,type="n",xlim=c(-2,12),ylim=c(-2,12))
#Plot on average 7 tree crowns of Weibull-distributed radius at random locations
n<-rpois(1,7)
circle(x=runif(n,0,10),
       y=runif(n,0,10),
       r=rweibull(n,6,2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ circle } 
% __ONLY ONE__ keyword per line
