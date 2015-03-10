\name{HD models}
\alias{HDnaslund}
\alias{HDmichailoff}
\alias{HDcurtis} 
\alias{HDmeyer}
\alias{HDpower}
\alias{HDnaslund2} 
\alias{HDnaslund3} 
\alias{HDnaslund4} 
\alias{HDmicment} 
\alias{HDmicment2} 
\alias{HDwykoff} 
\alias{HDprodan} 
\alias{HDlogistic} 
\alias{HDrichards} 
\alias{HDweibull} 
\alias{HDgomperz} 
\alias{HDsibbesen} 
\alias{HDkorf} 
\alias{HDratkowsky} 
\alias{HDhossfeldIV}
\alias{startHDnaslund}
\alias{startHDmichailoff}
\alias{startHDmeyer}
\alias{startHDpower}
\alias{startHDnaslund2} 
\alias{startHDnaslund3} 
\alias{startHDnaslund4} 
\alias{startHDmicment} 
\alias{startHDmicment2} 
\alias{startHDwykoff} 
\alias{startHDcurtis} 
\alias{startHDprodan} 
\alias{startHDlogistic} 
\alias{startHDrichards} 
\alias{startHDweibull} 
\alias{startHDgomperz} 
\alias{startHDsibbesen} 
\alias{startHDkorf} 
\alias{startHDratkowsky} 
\alias{startHDhossfeldIV}
\alias{startHDodels}
\alias{HDmodels}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Available 2- and 3- parameter H-D model functions to be used by function \code{fithd}.
}
\description{
Nonlinear functions for modeling tree height on diameter. Usually called using \code{fithd}. 
}
\usage{HDnaslund(d, a, b, bh=1.3) 
HDcurtis(d, a, b, bh=1.3) 
HDmichailoff(d, a, b, bh=1.3) 
HDmeyer(d, a, b, bh=1.3) 
HDpower(d, a, b, bh=1.3)
HDnaslund2(d, a, b, bh=1.3)
HDnaslund3(d, a, b, bh=1.3)
HDnaslund4(d, a, b, bh=1.3)
HDmicment(d, a, b, bh=1.3) 
HDmicment2(d, a, b, bh=1.3) 
HDwykoff(d, a, b, bh=1.3) 


HDprodan(d, a, b, c, bh=1.3) 
HDlogistic(d, a, b, c, bh=1.3) 
HDrichards(d, a, b, c, bh=1.3) 
HDweibull(d, a, b, c, bh=1.3) 
HDgomperz(d, a, b, c, bh=1.3) 
HDsibbesen(d, a, b, c, bh=1.3) 
HDkorf(d, a, b, c, bh=1.3) 
HDratkowsky(d, a, b, c, bh=1.3) 
HDhossfeldIV(d, a, b, c, bh=1.3)

startHDnaslund(d, h, bh=1.3) 
startHDcurtis(d, h, bh=1.3) 
startHDmichailoff(d, h, bh=1.3) 
startHDmeyer(d, h, bh=1.3) 
startHDpower(d, h, bh=1.3)
startHDnaslund2(d, h, bh=1.3) 
startHDnaslund3(d, h, bh=1.3) 
startHDnaslund4(d, h, bh=1.3) 
startHDmicment(d, h, bh=1.3) 
startHDmicment2(d, h, bh=1.3) 
startHDwykoff(d, h, bh=1.3) 

startHDprodan(d, h, bh=1.3) 
startHDlogistic(d, h, bh=1.3) 
startHDrichards(d, h, bh=1.3, b=0.04) 
startHDweibull(d, h, bh=1.3) 
startHDgomperz(d, h, bh=1.3) 
startHDsibbesen(d, h, bh=1.3, a=0.5) 
startHDkorf(d, h, bh=1.3) 
startHDratkowsky(d, h, bh=1.3, c=5) 
startHDhossfeldIV(d, h, bh=1.3, c=5)
}

\arguments{
  \item{d}{
  A vector of tree diameters, usually in cm
  }
  \item{h}{
  A vector of tree heights, usually in m. The observed heights should be always above or equal to \code{bh}. 
  }
  \item{a, b, c}{
  Parameters a, b (and c for 3- parameter functions) of the applied function. See details for expressions of different functions.
  }
  \item{bh}{
  The applied height for the measurement of tree diameter (so called breast height). Of the same unit as \code{h}.  
  }
  }

\details{
The available 2- parameter functions are
\itemize{
\item Naslund:          \eqn{h(d) = bh + \frac{d^2}{(a + bd)^2}}{h(d) = bh + d^2/(a + b d)^2}
\item Curtis:            \eqn{h(d) = bh + a \left(\frac{d}{1 + d}\right)^b}{h(d) = bh + a (d/(1 + d))^b}
\item Michailoff:        \eqn{h(d) = bh + a e^{-b d^{-1}}}{bh + a e^(-b d^(-1))}
\item Meyer:             \eqn{h(d) = bh + a (1-e^{-b d})}{h(d) = bh + a (1-exp(-b d))}  
\item Power:             \eqn{h(d) = bh + a d^b} 
\item Naslund2:          \eqn{h(d) = bh + \frac{d^2}{\left(a + e^b d\right)^2}}{h(d) = bh + d^2/(a + exp(b) d)^2}
\item Naslund3:          \eqn{h(d) = bh + \frac{d^2}{(e^a + b d)^2}}{h(d) = bh + d^2/(exp(a) + b d)^2}
\item Naslund4:          \eqn{h(d) = bh + \frac{d^2}{(e^a + e^b d)^2}}{h(d) = bh + d^2/(exp(a) + exp(b) * d)^2}
\item Michaelis-Menten:  \eqn{h(d) = bh + \frac{a d}{b + d}}{h(d) = bh + a d/(b + d)}
\item Michaelis-Menten2: \eqn{h(d) = bh + \frac{d}{a + b * d}}{h(d) = bh + d/(a + b * d)}
\item Wykoff:            \eqn{h(d) = bh + \exp\left(a + \frac{b}{d + 1}\right)}{h(d) = bh + exp(a + b/(d + 1))}
}
 
The available 3- parameter functions are
\itemize{
\item Prodan:           \eqn{h(d) = bh + \frac{d^2}{a + bd + c d^2}}{h(d) = bh + d^2/(a + bd + c d^2) }
\item Logistic:         \eqn{h(d) = bh + \frac{a}{1 + b e^{-c d}}}{h(d) = bh + a/(1 + b exp(-c d)) }
\item Chapman-Richards: \eqn{h(d) = bh + a (1 - e^{-bd})^c }{h(d) = bh + a (1 - exp(-bd))^c }  
\item Weibull:          \eqn{h(d) = bh + a (1 - e^{-b d^c}) }{h(d) = bh + a (1 - exp(-b d^c)) }
\item Gomperz:          \eqn{h(d) = bh + a \exp(-b \exp(-c d)) }{h(d) = bh + a exp(-b exp(-c d)) }
\item Sibbesen:         \eqn{h(d) = bh + a d^{b d^{-c}}}{h(d) = bh + a d^(b d^(-c))}
\item Korf:             \eqn{h(d) = bh + a \exp(-b d^{-c}) }{h(d) = bh + a exp(-b d^(-c)) }
\item Ratkowsky:        \eqn{h(d) = bh + a \exp\left(\frac{-b}{d + c}\right)}{h(d) = bh + a exp(-b/(d + c))}  
\item Hossfeld IV:      \eqn{h(d) = bh + \frac{a}{1 + \frac{1}{bd^c}}}{h(d) = bh + a/(1 + 1/(b*d^c))}
}

For each model, two functions are provided: one computing the value of the H-D model for given diameters using given values of parameters a, b (and c),
and another returning the initial guesses of a, b (and c) for given h-d data.  

The initial guesses are in most cases computed by fitting a linearized version of the model into the provided h-d data using \code{lm}. 
For some 3- parameter versions, 
no straightforward linearization is possible and one of the parameters is set to a fixed sensible constant. 
Those values can be seen as additional arguments in the corresponding startHD - functions. 
Details can be seen directly from the function definitions.

The user can define her own functions to be used with \code{fithd}. 
The case-sensitive naming of the functions should follow exactly the naming convention 
shown above. In addition, the names of the of arguments, as well as their order, 
should be the same as in the functions above. 

The models are named according to references in 
\itemize{
\item Zeide, B. 1993. Analysis of growth equations. Forest Science 39(3):594-616.
\item Huang, S., Titus, S.J., and Wiens, D.P. 1992. Comparison of nonlinear height-diameter functiond for major Alberta tree species. Can J. For. Res. 22: 1297-1304. 
}

Suggestions on naming and references on the functions are welcome. 
}

\value{
For functions HDxxx, a vector of tree heights corresponding diameters \code{d} is returned. 
For functions startHDxxx, a named vector of initial estimates of a, b and (c).
}

\author{
Lauri Mehtatalo lauri.mehtatalo@uef.fi
}
\examples{
data(spati)
theta<-startHDnaslund(spati$d,spati$h)
plot(spati$d,spati$h)
d<-seq(0,50)
lines(d,HDnaslund(d,a=theta[1],b=theta[2]),col="red",lwd=5)
}

\references{
Mehtatalo, L., Gregoire, T.G., and de Miguel, S. Modeling Height-diameter curves for height prediction. Manuscript.

Mehtatalo, L. 2004. A longitudinal height-diameter model for Norway spruce in Finland. Canadian Journal of Forest Research 34(1): 131-140.

Mehtatalo, L. 2005. Height-diameter models for Scots pine and birch in Finland. Silva Fennica 39(1): 55-66.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
