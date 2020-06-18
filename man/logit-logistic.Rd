\name{logit-logistic}
\alias{dll}
\alias{pll}
\alias{qll} 
\alias{rll}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Four-parameter Logit-logistic Distribution
}
\description{
Density, distribution function, quantile function and random generation for the four-parameter logit-logistic distribution. 
}
\usage{dll(x, mu, sigma, xi=0, lambda=1, log = FALSE)
pll(q, mu, sigma, xi=0, lambda=1, lower.tail=TRUE, log.p=FALSE)
qll(p, mu, sigma, xi=0, lambda=1, lower.tail=TRUE, log.p=FALSE)
rll(n, mu, sigma, xi=0, lambda=1)
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
  \item{mu,sigma,xi,lambda}{
  parameters of the distribution, \code{xi} (minimum) defaults to 0 and \code{lambda} (max-min) to 1. 
  }
  \item{log, log.p}{
  logical; if TRUE, probabilities p are given as log(p)  
  }
  \item{lower.tail}{
  logical; if TRUE (default), probabilities are \eqn{P[X \le x]} otherwise \eqn{P[X > x]}.
  }
  }

\details{
The logit-logistic cdf and pdf are
\deqn{F(d|\xi,\lambda,\mu,\sigma)  = 
 \frac{1}{1+e^{(\frac{\mu}{\sigma})}
  (\frac{d-\xi}{\xi+\lambda-d})^{-\frac{1}{\sigma}}}}
  
\deqn{f(d|\xi,\lambda,\mu,\sigma) =
  \frac{\lambda}{\sigma}\frac{1}{(d-\xi)(\xi+\lambda-d)}
  \frac{1}{e^{-\frac{\mu}{\sigma}}(\frac{d-\xi}{\xi+\lambda-d})^{\frac{1}{\sigma}}+e^{\frac{\mu}{\sigma}}(\frac{d-\xi}{\xi+\lambda-d})^{-\frac{1}{\sigma}}+2}}
Parameter \eqn{\xi} is the minimum, \eqn{\lambda>0} the width of range (max-min), \eqn{\mu} controls the skewness and \eqn{\sigma} the curtosis.  
}

\value{
\code{dll} gives the density, \code{pll} gives the distribution function, \code{qll} gives the quantile function, and \code{rll} generates random deviates.

Invalid arguments will result in return value \code{NaN}.

The length of the result is determined by \code{n} for \code{rll}, and is the maximum of the lengths of the numerical arguments for the other functions.

The numerical arguments other than \code{n} are recycled to the length of the result. Only the first elements of the logical arguments are used. 
}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
}
\examples{
data(spati)
d<-spati$d[spati$plot==22]
hist(d,freq=FALSE)
d0<-seq(0,60,0.1)
lines(d0,dll(d0,0.630,0.573,3.561,35.2))
}
\references{
Mingliang Wang and Keith Rennolls, 2005. Tree diameter distribution modelling: introducing the logit-logistic distribution. Canadian Journal of Forest Research, 35(6): 1305-1313, \doi{10.1139/x05-057}.

Mehtatalo, Lauri and Lappi, Juha 2020. Biometry for Forestry and Environmental Data: with examples in R. 
New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ forest biometrics}
