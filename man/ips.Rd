\name{ips}
\alias{ips}
\docType{data}
\title{
Wood-decaying fungi carried by bark beetle individuals and their mites.
}
\description{
The total number of fungal species (ophistomatoid and non-ophistomatoid fungi are coded separately) 
associated with \emph{Ips typographus} bark beetle individuals and their mites.
}
\usage{data(ips)}
\format{
     A data frame with 298 observations (bark betle individuals) on the following 5 variables.
    \describe{
    \item{\code{Fungi}}{The total number of fungal species associated with the individual bark beetle.}
    \item{\code{Ophi}}{The number of ophistomatoid fungal species.}
    \item{\code{Other}}{The number of non-ophistomatoid fungal species. The three first variables are related through \eqn{Other+Ophi=Fungi}.}
    \item{\code{Season}}{Categorical time of data collection with three levels: spring, summer or fall. The default is spring.}
    \item{\code{Mites}}{The number of mites found in the bark beetle.}
    }
}
\details{
The ophiostomatoid fungal families \emph{Microascales} and
\emph{Ophiostomatales} are common associates of bark beetle \emph{Ips
typographus}, which they use to spread within the wood material. The number of
fungal species in these families is high, and a certain beetle individual can 
carry several fungal species with it. The bark beetles may have
mites attached to them, and it may be possible that some
fungal species are associated to the beetles only through the mites. 

The dataset includes measurements of 289 bark beetle
individuals from a storm-felled Norway spruce forest in eastern Finland. For
each individual, the number of attached mites was determined using a
microscope. In addition the number of fungal species per
bark beetle was determined genetically. However, it was not possible to determine
whether the fungi were associated with the mites or the bark beetle itself. The
observations were collected at three different seasons: spring, summer and
fall of the same year, approximately 100 individuals in each season. The data are used to analyze the
effects of season and number of mites on the number of fungal species per bark
beetle. 
}

\references{Linnakoski, R., Mahilainen, S., Harrington, A., Vanhanen, H., Eriksson, M., Mehtatalo, L., 
Pappinen, A., Wingfield, M.J. 2016. The seasonal succession of fungi associated with Ips typographus 
beetles and their phoretic mites in an outbreak region of Finland. PLOS ONE. 10.1371/journal.pone.0155622. 
}

\examples{
data(ips)

ips$Mites2<-ips$Mites-mean(ips$Mites)

mod1<-glm(Fungi~Season+Mites,family=poisson,data=ips)
}
\keyword{datasets}
