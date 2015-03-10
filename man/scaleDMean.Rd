\name{scaleDMean}
\alias{scaleDMean}
\alias{scaleDGMean}
\alias{scaleDMed}
\alias{scaleDGMed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Weibull scale parameter for the given mean/median diameter and shape parameter. 
}
\description{
The function finds such scale parameter of the Weibull distribution that 
yields the given mean/median diameter. Function scaleDMean is used for  
arithmetic mean, scaleDGMean for the mean of basal-area weighted distribution, 
scaleDMed for median and scaleDGMed for the median of the basal-area weighted 
diameter distribution. 

The functions are used in the recovery of Weibull parameters using function 
\code{\link{recweib}}.
}
\usage{
scaleDMean(D,shape)
scaleDGMean(D,shape)
scaleDMed(D,shape)
scaleDGMed(D,shape)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D }{The diameter}
\item{shape}{The Weibull shape parameter}
}

\value{
\item{scale}{The value of the Weibull scale parameter.}
}
\references{
Siipilehto, J. and Mehtatalo, L. 2013. Parameter recovery vs. parameter prediction for the
Weibull distribution validated for Scots pine stands in
Finland. Silva Fennica 47(4), article id 1057.
}

\author{
Lauri Mehtatalo and Jouni Siipilehto
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{recweib}}
}
\examples{
scaleDMean(15,3)
scaleDGMean(15,3)
scaleDMed(15,3)
scaleDGMed(15,3)
}

\keyword{ forest biometrics }
