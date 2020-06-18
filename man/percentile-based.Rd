\name{percentile-based}
\alias{dPercbas}
\alias{pPercbas}
\alias{qPercbas} 
\alias{rPercbas}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Percentile-based Distribution
}
\description{
Density, distribution function, quantile function and random generation for the percentile-based distribution. 
}
\usage{dPercbas(x, xi, F)
pPercbas(q, xi, F)
qPercbas(p, xi, F)
rPercbas(n, xi, F)
}

\arguments{
  \item{x, q}{
  vector of quantiles
  }
  \item{p}{
  vector of probabilitiies
  }
  \item{n}{
  number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
  }
  \item{xi}{
  Strictly increasing vector of percentiles corresponding to the cumulative probabilities given in \code{F}. Of same length as \code{F}. 
  }
  \item{F}{
  a \code{k}-length strictly increasing  vector of cumulative probabilities, with \code{xi[1]=0} and \code{xi[k]=1}. 
  }
  }

\details{
The percentile-based distribution is defined by the quantiles \code{xi} that correspond to the cumulative probabilities given in \code{F}. 
The continuous distribution is obtained by linear interpolation of the cdf. 
}

\value{
\code{dll} gives the density, \code{pll} gives the distribution function, \code{qll} gives the quantile function, and \code{rll} generates random deviates.

The length of the result is determined by \code{n} for \code{rPercbas}, and by the length of \code{x}, \code{q} and \code{p} for the other functions.
}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
}
\examples{
d0<-seq(0,30,0.01)
d<-c(5,7,10,11,11.7,13,15,19,22,24,25,28.5)
plot(d0,pPercbas(d0,d),type="l")
hist(rPercbas(1000,d),breaks=seq(0,30,1),freq=FALSE,ylim=c(0,0.15))
lines(d0,dPercbas(d0,d),col="red")
}
\references{
Borders B. E., Souter R. A., Bailey. R. L., and Ware, K. D. 1987. Percentile-based distributions characterize forest stand tables. Forest Science 33(2): 570-576.

Mehtatalo, L. 2005. Localizing a predicted diameter distribution using sample information. Forest Science 51(4): 292--302. 

Mehtatalo, Lauri and Lappi, Juha 2020a. Biometry for Forestry and Environmental Data: with examples in R. New York: Chapman and Hall/CRC. 426 p. 
\doi{10.1201/9780429173462}

Mehtatalo, Lauri and Lappi, Juha 2020b. Biometry for Forestry and Environmental Data: with examples in R. Full Versions of The Web Examples. Available at \url{http://www.biombook.org}. 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ forest biometrics}
