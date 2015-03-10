\name{circle}
\alias{circle}
\title{Plot circles of a specified radius}
\description{Adds circles of radii r at coordinates specified by x and y onto an existing plot.}
\usage{
circle(x, y, r, col = "black", lty = "solid", lwd = 1, grayfill = FALSE)
}
\arguments{
  \item{x,y}{the x- and y- coordinates of the midpoints}
  \item{r}{the radii of the circles}
  \item{col, lty, lwd}{the draving color, line type and line width of the perimeter line}
  \item{grayfill}{Logical- whether gray shading will be used for the circles or not
}

}
\value{
This function is used for its side effects.
}

\author{
Lauri Meht\"atalo
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
