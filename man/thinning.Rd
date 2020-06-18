\name{thinning}
\alias{thinning}
\docType{data}
\title{
Effect of thinning on individual tree growth
}
\description{
A time series of estimated effect of thinning on the annual basal area growth of a 
Scots pine tree. 
}
\usage{data(thinning)}
\format{
  A data frame with 23 observations on the following 3 variables.
  \describe{
    \item{\code{TreeID}}{Tree ID, 3_3 for all observations in this data.}
    \item{\code{Year}}{Calendar year (1983-2005).}
    \item{\code{ThEff}}{Estimated effect of thinning on the annual basal area growth in mm^2.}
    }
}
\details{
The thinning took place between years 1986 and 1987. For details about the original measurements, see the documentation of data set \link{patti}, \link{afterthin}. For details about the procedure used 
in extracting the thinning effect, see Example 6.6 in Mehtatalo and Lappi 2020b and for nonlinear modeling of this data, see Chapter 7 of Mehtatalo and Lappi 2020a.}

\references{
Mehtatalo, L., Peltola, H., Kilpelainen, A. and Ikonen, V.-P. 2014. The response of basal area growth 
of Scots pine to thinning: A longitudinal analysis of tree-specific series using a nonlinear 
mixed-effects model. Forest Science 60 (4): pp. 636-644. \doi{10.5849/forsci.13-059}.

Mehtatalo, Lauri and Lappi, Juha 2020a. Biometry for Forestry and Environmental Data: with examples in R. 
New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}

Mehtatalo, Lauri and Lappi, Juha 2020b. Biometry for Forestry and Environmental Data: with examples in R. 
Full Versions of The Web Examples. Available at \url{http://www.biombook.org}. 
}

\seealso{
  \code{\link{patti}}, \code{\link{afterthin}}, \code{\link{thefdata}}.
}
\examples{
data(thinning)
plot(thinning$Year,thinning$ThEff,type="l")
}
\keyword{datasets}
