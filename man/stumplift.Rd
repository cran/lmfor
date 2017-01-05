\name{stumplift}
\alias{stumplift}
\docType{data}
\title{
Productivity of stump lifting machines.
}
\description{
The productivity of stump lifting machines on three Norway Spruce (Picea Abies) clearcut areas (sites). 
Stumps are lifted for use as bioenergy. The data were collected from three sites in 
Central Finland.  
}
\usage{data(stumplift)}
\format{
  A data frame with 485 observations on the following 5 variables.
  \describe{
    \item{\code{Stump}}{A unique stump id based on the order of processing. The successive numbers 
    are usually close to each other in the clearcut area, but nearby trees do not necessarily 
    have small difference in stump id.}
    \item{\code{Machine}}{The machine/clearcut/dirver combination. A factor with three levels.}
    \item{\code{Diameter}}{Stump diameter, cm.}
    \item{\code{Time}}{Processing time, seconds.}
    \item{\code{Productivity}}{Productivity, \eqn{m^3}/effective working hour}
  }
}
\details{
Each site was operated with different machine and driver so that the 
effect of site, machine and driver cannot be separated. The volume of each stump was 
estimated using the function of Laitila (2008), based on the stump diameter. 
A work system study was conducted to measure the processing time (seconds) and 
productivity (\eqn{m^3}/hour) for each stump. 
}

\references{
Teijo Palander, Kalle Karha, Lauri Mehtatalo 2016. 
Applying polynomial regression modeling to productivity analysis of 
sustainable stump harvesting. Scandinavian Journal of Forest Reseach (Accepted). 
doi:10.1080/02827581.2016.1238957


Teijo Palander, Janne Smolander, Kalle Karha, 2015. 
Work system study of three stump-lifting devices in Finland. 
Scandinavian Journal of Forest Research 30(6) 558-567, doi:10.1080/02827581.2015.1027731
}

\examples{
data(stumplift)
library(nlme)

modConstPow<-gls(Productivity~Machine+Machine*I((Diameter-70)^2),
                 data=stumplift,
                 weights=varPower(),
                 corr=corAR1(form=~Stump|Machine))
           
}
\keyword{datasets}
