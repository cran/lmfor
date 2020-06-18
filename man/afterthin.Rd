\name{afterthin}
\alias{afterthin}
\docType{data}
\title{
Increment core data of Scots pine trees
}
\description{
Post-thinning growth ring measurements of 88 trees of a long-term thinning 
experiment on a naturally regenerated
Scots pine stand in Eastern Finland. 
}
\usage{data(afterthin)}
\format{
  A data frame with 1319 observations on the following 7 variables.
  \describe{
    \item{\code{Plot}}{Sample plot id, a factor with 10 levels.}
    \item{\code{Tree}}{Tree id, a factor with 55 levels (same tree id may occur on different plots!).}
    \item{\code{Year}}{Calendar year of the ring.}
    \item{\code{SDAfterThin}}{Stand density (trees per ha) of the sample plot.}
    \item{\code{SDClass}}{Thinning treatment, factor with 4 levels (1=Control, 2=Light, 3=Moderate, 4=Heavy).}
    \item{\code{CA}}{Current tree age in years.}
    \item{\code{RBA}}{Ring Basal area, \eqn{mm^2} }
  }
}
\details{
Long-term thinning experiment on a naturally regenerated
Scots pine stand in Eastern Finland. The experiment consists of 10
sample plots, in four different classes according
to the post-thinning stand density. The plots were thinned in winter 1986-1987. 
In winter 2006 -2007, 10 trees were felled from each plot.
A radial 5mm by 5mm segment from pith to bark was cut from each tree at height 1.3 meter height.
Ring widths from pith to bark were analyzed for
each sample, using an ITRAX X-ray microdensitometer an post-processed to
create ring widths from pith to bark were determined for
each disc. The ring widths were further transformed to ring basal areas by assuming circular, 
growth rings. For 12 trees, ring widths could not be extracted. The data includes 
ring widths for a total of 88 trees between years 1991-2005. The original data is available in data set 
\link{patti}.}

\references{
Mehtatalo, L., Peltola, H., Kilpelainen, A. and Ikonen, V.-P. 2014. The response of basal area growth 
of Scots pine to thinning: A longitudinal analysis of tree-specific series using a nonlinear 
mixed-effects model. Forest Science 60 (4): pp. 636-644. DOI: \doi{10.5849/forsci.13-059}.

Mehtatalo, Lauri and Lappi, Juha 2020. Biometry for Forestry and Environmental Data: with examples in R. 
New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}
}

\seealso{
  \code{\link{patti}}, \code{\link{thefdata}}, \code{\link{thinning}}.
}
\examples{
data(afterthin)
par(mfcol=c(2,1),cex=0.7,mai=c(0.8,0.8,0.5,0.1))
linesplot(afterthin$CA,
		  afterthin$RBA,
		  group=afterthin$Plot:afterthin$Tree,
		  col.lin=as.numeric(afterthin$SDClass),cex=0,
		  xlab="Tree age",
		  ylab=expression("Ring basal area, "*mm^2))
		  
linesplot(afterthin$Year,
		  afterthin$RBA,
		  group=afterthin$Plot:afterthin$Tree,
		  col.lin=as.numeric(afterthin$SDClass),cex=0,
		  xlab="Year",
		  ylab=expression("Ring basal area, "*mm^2))
}
\keyword{datasets}
