\name{BrkRes}
\alias{BrkRes}
\docType{data}
\title{
Breaking resistance (=bending strength) of birch wood samples
}
\description{
Measurements of the breaking resistance of
wood samples from a total of 118 downy birch (\emph{Betula pubescens}) trees. 
}
\usage{data(BrkRes)}
\format{
  A data frame with 274 observations on the following 7 variables.
  \describe{
    \item{\code{Tree}}{Tree id (numeric)}
    \item{\code{Resistance}}{Breaking resistance, MPa}
    \item{\code{Density}}{Wood density (as air-dry in 12-15\% moisture) , g/cm^3}
    \item{\code{FibreLength}}{Fibre length, mm}
    \item{\code{RingClass}}{Categorical with three levels indicating the position within the stem. 1=near the pitch (inside), 2=middle, 3=near the bark (outer).}
    \item{\code{SeedOrigin}}{Binary variable about the origin of the tree: 1=from seed, 0=sprouted. }
    \item{\code{Site}}{Categorical site class with four levels.}
  }
}
\details{
A
total of 1 - 4 wood samples per tree were collected.
The samples were classified to three ring classes according to the distance of the sample
from the pith. Each sample was measured
destructively for breaking resistance in the laboratory. The measured variables 
also included wood density. The data has been collected by Hanna Joronen, Katri 
Luostarinen and  Veikko Mottonen.
}

\references{
Mehtatalo, Lauri and Lappi, Juha 2020. Biometry for Forestry and Environmental Data: with examples in R. 
New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}

Joronen, H. 2020. Taivutusmurtolujuuteen ja kimmokertoimeen vaikuttavat tekijat hieskoivulla 
nuorpuussa ja aikuispuussa seka niiden mallintaminen lineaarisella sekamallilla. Master's thesis, University of Eastern Finland. 
}

\examples{
data(BrkRes)

brmod1 <- lme(Resistance~RingClass+Density, 
          random=~RingClass-1|Tree, 
          data=BrkRes)
}
\keyword{datasets}
