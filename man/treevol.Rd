\name{treevol}
\alias{treevol}
\docType{data}
\title{
Individual tree volume modeling data for 8508 pine, spruce and birch trees in Finland
}
\description{
Measurements of tree volumes from 8508 individual Scots pine (4066 trees), 
Norway spruce (2966 trees) and birch (1476 trees) trees in Finland. 
A total of 5053 trees were measured in 1970's by climbing to the standing trees, 
1534 trees by falling the trees in 1990's, 
and 1921 trees by using terrestrial laser scanning in 2010's. The trees have been 
selected from among the sample trees of national forest inventories. 
In addition to tree volumes, diameters, heights and upper diameters, the data set includes 
variables describing the grouping of the trees to stands and plots and information about the site characteristics. 
}
\usage{data(treevol)}
\format{
  A data frame with 8508 observations on the following 11 variables.
  \describe{
    \item{\code{tree_id}}{Tree id (numeric)}
    \item{\code{stand}}{Stand id (numeric)}
    \item{\code{plot}}{Plot id (nested within stand)}
    \item{\code{v}}{Tree volume, dm^3}
    \item{\code{dbh}}{Tree diameter at 1.3m height above the stump height, cm}
    \item{\code{h}}{Tree height, m}
    \item{\code{d6}}{Tree diameter at 6m height above the stump height, cm}
    \item{\code{dataset}}{Categorical indicator for the dataset, (see the general description).}
    \item{\code{temp_sum}}{Numeric temperature sum of the site. A high value indicates warm site. }
    \item{\code{species}}{Categorical with three levels, tree species.}
    \item{\code{soil}}{Categorical with two levels: mineral soil or peatland.}
    \item{\code{region}}{Categorical with four levels: region.}
  }
}
\details{
The data combines three datasets: 1) The climbed data including 2326 
Scots pine (Pinus sylvestris L.), 1864 Norway spruce (Picea abies L. Karst.) 
and 863 silver (Betula pendula Roth) and downy (Betula pubescens Ehrh.) 
birches collected 1968-1972, and 2) the felled data with 797 pines, 479 spruces 
and 258 birches collected 1988-2001, and the scanned data with 943 Scots pine, 
623 Norway spruce and 355 downy or silver birch collected 2017-2018. 
For the two older dat sets, diameter measurements are based on two perdendicular measurements 
of the stem at regularly spaced height. In the scanned data, the diemater measurements 
are based on TLS but heights are measured manually. 
For more details concerning the climbed data, see Laasasenaho (1982), 
and for the scanned data see Pitkanen et al.  (2021) and for the 
felled data see (Korhonen & Maltamo 1990). 
In the climbed and felled data, the volumes of the trees were predicted using a natural 
cubic spline (splinefun, R Core Team 2021) to interpolate the taper curve between 
measured diameters at given relative heights. The volumes were obtained as a slolid of revolution based on the 
interpolated taper curve.
In scanned data,  a smoothing spline (smooth.spline, R core team 2021) 
function was used instead of an interpolating spline.
The pooled data has been used in Kangas et al. (XXXX). 
}

\references{
Kangas A., Pitkanen T., Mehtatalo L., Heikkinen J. (xxxx)
Mixed linear and non-linear tree volume models with regionally varying parameters

Korhonen, K. and Maltamo, M. 1990. Mannyn maanpaallisten osien kuivamassat Etela-Suomessa. 
Metsantutkimuslaitoksen tiedonantoja 371. 

Laasasenaho, J. (1982). Taper curve and volume functions for pine, spruce and birch. 
Communicationes Instituti Forestalia Fennica 108: 1-74.

Pitkanen, T.P., Raumonen, P., Liang, X., Lehtomaki, M. and Kangas, A.S. 2020. 
Improving TLS-based stem volume estimates by field measurements. 
Computers and Electronics in Agriculture and Forestry 180, 105882.
}

\examples{
\dontrun{
data(treevol)
        
treevol$formfactor<-treevol$v/volvff(treevol$dbh,treevol$h,logita=100,lambda=log(0.2))
treevol$logitff<-log((treevol$formfactor)/(1-(treevol$formfactor)))

ptrees<-treevol[treevol$species=="pine",]
mod.init<-lm(logitff~I(1/h)+h+dbh+I(h*dbh)+I(1/(h*dbh))+
                     dataset+dataset:dbh+dataset:h,data=ptrees)

mod<-nlme(v~volvff(dbh,h,logita=logita,lambda=lambda),
         fixed=list(logita~I(1/h)+h+dbh+I(h*dbh)+I(1/(h*dbh))+dataset+
                           dataset:dbh+dataset:h+soil+temp_sum,
                           lambda~1),
         random=logita~1|stand/plot,
         start=c(coef(mod.init),rep(0,2),log(0.2)),
         data=ptrees,
         weights=varComb(varIdent(form=~1 |dataset),varPower()),
         method="ML",
         control=list(msVerbose=TRUE),
         verbose=TRUE)
}                    
}
\keyword{dataset}
\keyword{tree volume}
