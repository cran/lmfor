\name{spati}
\alias{spati}
\title{Raw sample plot data of Scots pine in Ilomantsi, Finland.}
\description{A dataset of Scots pine growth. 
The trees were collected on 56 fixed-area sample plots. The data includes no remeasurements. 
The growth data are based on measurements of increment borer chips.}
\usage{data(spati)}
\format{
  A data frame with 9913 observations on the following 3 variables.
  \describe{
    \item{\code{plot}}{A unique sample plot id.}
    \item{\code{X}}{x-coordinate of the plot}
    \item{\code{Y}}{Y-coordinate of the plot}
    \item{\code{N}}{Stand density, trees per ha}
    \item{\code{G}}{Basal area, \eqn{m^2/ha}}
    \item{\code{V}}{Plot volume, \eqn{m^3/ha}}
    \item{\code{Dg}}{Basal-area weighted mean diameter, cm} 
    \item{\code{Hg}}{Height of basal area median diameter tree, m}
    \item{\code{Tg}}{Age of basal area median tree, yr}
    \item{\code{Hdom}}{Dominant height, m}
    \item{\code{maos}}{percentage of Scots pines of the total volume}
    \item{\code{kuos}}{percentage of Norway spruces of the total volume}
    \item{\code{kanro}}{A unique sample plot id (same as plot).}
    \item{\code{puunro}}{Tree id within plot.}
    \item{\code{pl}}{tree species. 1=Scots Pine}
    \item{\code{xk}}{x- coordinates of trees within plot}
    \item{\code{yk}}{y- coordinates of trees within plot}
    \item{\code{d}}{Tree diameter at breast height (1.3 meters above the ground) in cm.}
    \item{\code{h}}{Tree height, m.}
    \item{\code{t}}{Tree age, years}
    \item{\code{dk}}{Tree diameter at stump height, cm. there seems to be some unclear issues.}
    \item{\code{X2b}}{Double bark thickness, mm}
    \item{\code{id1}}{Tree diameter growth within the 5 year period prior to the measurement. Missing data coded as -1.}
    \item{\code{id2}}{Tree diameter growth within the period 6-10 years prior to the measurement. Missing data coded as -1.}
  }
}

\author{
The data were collected by Timo Pukkala.
}

\references{
Pukkala, T. 1989. Prediction of tree diameter and height in a Scots pine stand as a
function of the spatial pattern of trees. Silva Fennica 23(2): 83-99.
}

\keyword{ datasets }
