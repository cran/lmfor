\name{scaleDMean}
\alias{scaleDMean1}
\alias{scaleDGMean1}
\alias{scaleDMed1}
\alias{scaleDGMed1}
\alias{scaleDMean2}
\alias{scaleDGMean2}
\alias{scaleDMed2}
\alias{scaleDGMed2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Weibull scale parameter for the given mean/median diameter and shape parameter. 
}
\description{
The function finds such scale parameter of the Weibull distribution that 
yields the given mean/median diameter. Function scaleDMean is used for  
arithmetic mean, scaleDGMean for the mean of basal-area weighted distribution, 
scaleDMed for median and scaleDGMed for the median of the basal-area weighted 
diameter distribution. Functions with number 1 in the name use Weibull functions 
as the unweighted density and functions with value 2 in the name use Weibull function as the 
basal-area weighted density. 

The functions are used in the recovery of Weibull parameters using function 
\code{\link{recweib}}.
}
\usage{
scaleDMean1(D,shape)
scaleDGMean1(D,shape)
scaleDMed1(D,shape)
scaleDGMed1(D,shape)
scaleDMean2(D,shape)
scaleDGMean2(D,shape)
scaleDMed2(D,shape)
scaleDGMed2(D,shape)
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
Finland. Silva Fennica 47(4), article id 1057. \doi{10.14214/sf.1057}
}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}> and Jouni Siipilehto
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{recweib}}
}
\examples{
scaleDMean1(15,3)
scaleDGMean1(15,3)
scaleDMed2(15,3)
scaleDGMed2(15,3)
}

\keyword{ forest biometrics }
\keyword{ diameter distributions }
