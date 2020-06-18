\name{thefdata}
\alias{thefdata}
\docType{data}
\title{
Effect of thinning on individual tree growth for 62 trees.
}
\description{
Time series of estimated effect of thinning on the annual basal area growth of 62 
Scots pine trees from three different thinning intensities.  
}
\usage{data(thefdata)}
\format{
  A data frame with 1238 observations on the following 7 variables.
  \describe{
    \item{\code{Plot}}{Sample plot id, a factor with 10 levels.}
    \item{\code{Tree}}{Tree id, a factor with 55 levels (same tree id may occur on different plots!).}
    \item{\code{Year}}{Calendar year of the ring.}
    \item{\code{SDClass}}{Thinning treatment, factor with 4 levels (1=Control, 2=Light, 3=Moderate, 4=Heavy).}
    \item{\code{CA}}{Current tree age in years.}
    \item{\code{Diam1986}}{Tree diameter in 1986 (just before the thinning).}
    \item{\code{ThEf}}{The thinning effect in annual basal area growth, mm^2.}
  }
}
\details{
The data are based on a long-term thinning experiment on a naturally regenerated
Scots pine stand in Eastern Finland, see documentation of \link{afterthin} and \link{patti} for details. 
The experiment consists of 10
sample plots, in four different classes (including an unthinned control) according
to the post-thinning stand density. The plots were thinned in winter 1986-1987. 
The undisturbed growth for each tree was predicted using a linear mixed-effect 
model with crossed year and tree effect, which was fitted to a data set including the control 
treatment for all calendar years and other trees until the year of thinning. 
The thinning effect was computed as the difference of observed growth and the undisturbed predicted growth. 
For details about the procedure used 
in extracting the thinning effect, see Example 6.6 in Mehtatalo and Lappi 2020b and for nonlinear modeling of this data, see Chapter 7 of Mehtatalo and Lappi 2020a.
The data includes only the observations after the thinning year for the tree thinned treatments (i.e., control is excluded). 
}

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
  \code{\link{patti}}, \code{\link{afterthin}}, \code{\link{thinning}}.
}
\examples{
data(thefdata)
linesplot(thefdata$Year,thefdata$ThEf, thefdata$Tree,col.lin=thefdata$SDClass)
}
\keyword{datasets}
