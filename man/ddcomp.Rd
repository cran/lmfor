\name{ddcomp}
\alias{ddcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate the fit of a tree diameter distribution
}
\description{
A function to compare the fit of the observed tree diameter data (d) 
to a specified diameter distribution (density). 
}
\usage{
ddcomp(d,density="dweibull",power=0,limits=seq(0,100),limitsd=limits,plot=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{d}{numeric vector of observed diameters} 
\item{density}{either a valid name for a probability density function in R or
               a vector of diameter class densities for diameter classes whose 
               limits are given in vector limitsd}
\item{power}{the weight used in error index. Value 2 gives BA weight, 0 (default) the unweighted}
\item{limits}{the diameter class limits to compute the error index} 
\item{limitsd}{see the description of argument \code{density}} 
\item{plot}{logical. Should a graph be produced to illustrate the ecdf of d  and the cdf corresponding to density}
\item{...}{additional arguments passed to function specified by a character-type density. 
e.g. Weibull shape and scale of if density="dweibull"}
}

\details{
The comparison is  done for mean, variance and standard deviation and shape. 
The shape is compared by computing the sum of absolute differences 
(error index) in densities for the observed data and predicted density in diameter classes
specified by "limits". The error index has therefore a value between 0 (complete match) and 
2 (complete mismatch). The error index is computed for the predicted density as such (\code{ei1})
and to a rescaled and switched density, which has exactly same mean and variance as the 
given diameter data (\code{e12}). 

The error index is calculated as the sum of variable \eqn{(f_{obs}-f_{pred})x^{power}} over the diameter classes, where \eqn{x} is 
the midpoint of the diaemeter class and \eqn{f_{obs}-f_{pred}} is the difference in predicted and observed frequency. 
By default, \eqn{power=0}.
}
\value{
A list of components 
\item{mudif}{The difference in means}
\item{vardif}{The difference in variances}
\item{sddif}{The difference in standard deviations}
\item{ei1}{the error index for original predicted distribution (see details)}
\item{ei1}{the error index for scaled predicted distribution (see details)}
}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
}

\examples{
# Example 
# Observed diameters 
d<-c(18.8,24.2,18.7,13.0,18.9,22.4,17.6,22.0,18.8,22.9,
     16.7,13.7,20.6,15.1,31.8,17.2,19.6,16.8,19.3,27.4,
	 23.7,18.2,19.7,18.9,23.0,21.4,23.8,22.1,24.2,20.9)

# Weibull(5,20) distribution in 1 cm classes (class limits from 0,...,60)
f<-pweibull(1:60,5,20)-pweibull(0:59,5,20)

# compare using the classified true distribution (approximate)
ddcomp(d,density=f,limitsd=0:60,limits=0:100,plot=TRUE)

# compare b specifying a Weibull dsitribution (accurate)
ddcomp(d,density="dweibull",shape=5,scale=20,plot=TRUE)
}

\keyword{forest biometrics}
\keyword{diameter distributions}
