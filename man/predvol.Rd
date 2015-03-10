\name{predvol}
\alias{predvol}
\title{Individual tree volume functions for Finland}
\description{
Predict individual tree volumes using the functions of Laasasenaho(1982). 
The volume prediction can be based on tree diameter or tree diameter and height. 
The functions applying upper stem diameter have not (yet) been implemented.
}
\usage{
predvol(species,d,h=0,model=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{ The vector of tree species. 1:Pine, 2:Spruce, 3: Silver birch. 4: Downy birch.
                  Other codes than 1-4 are accepted but return NA as the volume prediction. 
}
  \item{d}{The vector of tree diameters at breast height (cm)
}
  \item{h}{
The vector of tree heights. Used only if model=2.
}
  \item{model}{
The model used. If model is 1, the prediction is based on tree diameter only. If model=2, then diameter and height are used.
}
}
\details{
Vectors species, dbh and height should be either scalars 
or vectors of the same length so that each element corresponds to 
one individual tree.
}
\value{
A vector of tree volumes (in liters).
}
\references{
Laasasenaho, Jouko 1982. Taper curve and volume functions for pine, spruce and birch. Comm. Inst. For. Fenn 108: 1-74.
}
\author{
Lauri Mehtatalo
}

\examples{
d<-c(15,18.3,29.3,22.4)
h<-c(13,18,22,19)
species<-c(1,1,1,3)
predvol(species,d,h,model=2)
predvol(species,d,model=1)
}