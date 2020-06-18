\name{linesplot}
\alias{linesplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A spaghetti plot of grouped data
}
\description{
Orders the observations by \code{x} and thereafter plots \code{y} on \code{x} and connects observations of the same group by lines. 
Useful, for example, to plot a longitudinal dataset.
}
\usage{
linesplot(x, y, group, xlab = "x", ylab = "y", 
          main = "", cex = 0.5, pch = 19, col = 1, col.lin = 1, 
          lw = FALSE, ylim = NULL, xlim = NULL, add = FALSE, lty = "solid", lwd=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{ Numerical vectors of the same length including the x and y variables.  
}
  \item{group}{
The variable specifying the group. Should be of the same length as vectors x and y.
}
  \item{xlab, ylab, 
        main, cex, 
        pch, col, 
        col.lin, xlim, 
        ylim, lty, lwd}{
Graphical parameters, see \code{\link{par}}
}
  \item{lw}{
Boolean. Whether a loess smoother to be added onto the plot.
}

  \item{add}{
Boolean. Whether to add to an existing plot or to open a new window. 
}
}
\details{
The observations within the group are connected at the increasing order of \code{x}. 
}
\value{
Used for its side effects.
}
\references{
Mehtatalo, Lauri and Lappi, Juha 2020. Biometry for Forestry and Environmental Data: with examples in R. 
New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}
}
\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
}

\examples{
D<-rep(seq(10,30),10)
H<-(20+rep(rnorm(10,0,0.5),each=21))*exp(-1.5*D^(-1.3))
plot<-rep(1:10,each=21)
linesplot(D,H,plot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ graphs} 
