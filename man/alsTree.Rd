\name{alsTree}
\alias{alsTree}
\docType{data}
\title{Individual tree characteristics and ALS data}
\description{
Field-measured and remotely sensed characteristics of 1510 individual Scots Pine trees from 56 sample plots in Kiihtelysvaara, Eastern Finland. 
}
\usage{data(alsTree)}
\format{
  A data frame with 1510 observations on the following 15 variables.
  \describe{
    \item{\code{plot}}{Sample plot id, integer}
    \item{\code{tree}}{Tree id, integer}
    \item{\code{DBH}}{Tree diameter at breast height, cm}
    \item{\code{H}}{Tree height, m}
    \item{\code{V}}{Tree volume, m^2)}
    \item{\code{HDB}}{Height of lowest dead branch, m}
    \item{\code{HCB}}{Crown base height, m}
    \item{\code{hmax}}{Maximum return height, m}
    \item{\code{h20}}{20th percentile of return heights within the tree crown, m}
    \item{\code{h30}}{30th percentile of return heights within the tree crown, m}
    \item{\code{h70}}{70th percentile of return heights within the tree crown, m}
    \item{\code{h80}}{80th percentile of return heights within the tree crown, m}
    \item{\code{a_hmean}}{Mean height of returns in the 250m^2 neighbourhood of the tree, m}
    \item{\code{a_veg}}{Proportion of returns from vegetation in the 250m^2 neighbourhood of the tree, m}
    \item{\code{a_h30}}{30th percentile of returns in the 250m^2 neighbourhood around the tree, m}
    \item{\code{a_h70}}{70th percentile of returns in the 250m^2 neighbourhood around the tree, m}
  }
}
\details{
Field measurements of tree diameter and height, height of dead branch and crown base height and tree location were taken from the trees in a field campaign. 
Volume was estimated based on diameter, height and upper stem diameter. 
In addition, the area was remotely sensed using airborne laser scanning. Detectable individual trees were delineated from the ALS point cloud and 
associated with the field measurements. From a large set of tree-specific ALS characteristics, the data includes those that were used in the final models of stand characteristics in Maltamo et al (2012).
}

\references{
Maltamo, M., Mehtatalo, L., Vauhkonen, J. and Packalen, P. 2012. Predicting and calibrating tree attributes 
by means of airborne laser scanning and field measurements. Canadian Journal of Forest Research 42: 1896-1907.
\doi{10.1139/x2012-134}

Mehtatalo, Lauri and Lappi, Juha 2020. Biometry for Forestry and Environmental Data: with examples in R. 
New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}.
}

\examples{
data(alsTree)
}

\keyword{datasets}
