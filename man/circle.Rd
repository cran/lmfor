\name{circle}
\alias{circle}
\title{Plot circles of a specified radius}
\description{Adds circles of radii r at coordinates specified by x and y onto an existing plot.}
\usage{
circle(x,y,r,border="black",lty="solid",lwd=1,fill=NULL)
}
\arguments{
  \item{x,y,r}{Vectors of the x- and y- coordinates of the midpoints and the associated radii. 
               Vectors x, y and r should be of the same length}
  \item{border, lty, lwd}{the draving color, line type and line width of the perimeter line. 
                          Use border=NA to omit the perimeter.}
  \item{fill}{The color used to fill the circles. fill=NULL does not fill at all.}

}
\value{
This function is used for its side effects on the graphical display.
}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
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
\keyword{ graphs} 
% __ONLY ONE__ keyword per line
