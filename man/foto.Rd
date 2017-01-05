\name{foto}
\alias{foto}
\docType{data}
\title{
CO2 exchange of transplanted Sphagnum fuscum moss in a chronosequence of mires.
}
\description{
The net carbon dioxide exchange of late successional moss species (Sphagnum fuscum) samples 
under seven levels 
of photosynthetic photon flux density in cronosequence of land uplift mires on the Finnish side 
of Bothnia Bay in Siikajoki, Finland. 
Moss samples were transplanted from the late succession site (Site 6) to all sites and photosynthetic activity 
was measured one year later for those samples which had survived.  
}
\usage{data(foto)}
\format{
  A data frame with 455 observations on the following 8 variables.
  \describe{
    \item{\code{Site}}{a factor with levels \code{1},...,6 from the earliest successional stage to the latest }
    \item{\code{Treatment}}{a factor with value \code{1} for samples with competitor removal treatment 
                            and \code{0} for untreated control.}
    \item{\code{sample}}{a factor with unique value for each of the 72 survived samples}
    \item{\code{moisture}}{moisture of the sample}
    \item{\code{PARtop}}{Photosyntetically active radiation (photon flux density  \eqn{\mu mol /m^2 /s^2})}
    \item{\code{WT}}{water table, cm}
    \item{\code{A}}{Net CO2 exchange, \eqn{\mu {mol} / g / h}, }
    \item{\code{subplot}}{a factor with unique value for each replicate}
  }
}
\details{
The number of transplanted 
replicates per site was 12, with two samples per replicate. One of the samples was treated with 
competing vegetation removal before transplanting whereas the other was left untreated. 
The 12 replicates per site were planted in locations with 2 to 3 different ground water table levels. 
A year after the transplanting, the photosynthetic activity (A) of the survived transplanted samples 
was recorded using seven artifically created light conditions ranging 
from complete darkness (PPFD=0) to extreme light conditions (PPFD=2000)  using an open, fully controlled flow-
through gas exchange fluorescence measurement system (GFS-3000; Walz, Effeltrich, Germany).
}

\references{
Laine, A.M., Ehonen, S., Juurola, E., Mehtatalo, L., and Tuittila, E-S. 2015. 
Performance of late succession species along a chronosequence: Environment does not exclude 
Sphagnum fuscum from the early stages of mire development. 
Journal of Vegetation Science 26(2): 291-301. DOI: 10.1111/jvs.12231
}

\examples{
data(foto)

LightResp<-function(PPFD,alpha=0.1,Pmax=10,A0=0) {
	       A0+Pmax*PPFD/(alpha+PPFD)
           }
               
library(nlme)

model5<-nlme(A~LightResp(PARtop,alpha,Pmax,A0),
		fixed=list(alpha~Site+Treatment+moisture,Pmax~Site+Treatment,A0~Site),
		random=list(sample=Pmax+alpha~1),
		data=foto,
		start=c(c(80,0,0,0,0,0,0,0),c(100,0,0,0,0,0,0),c(-20,0,0,0,0,0)),
		verbose=TRUE)
}
\keyword{datasets}
