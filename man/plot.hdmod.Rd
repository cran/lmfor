\name{plot.hdmod}
\alias{plot.hdmod}
%\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic plot a Height-Diameter model residuals
}
\description{
Plotting method for class \code{hdmod}
}
\usage{
\method{plot}{hdmod}(x, col.point = "blue", highlight = "red", standd = TRUE, 
           cex=1, corD=FALSE, ask=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A H-D model model fitted by \code{\link[lmfor]{fithd}}.
  }
  \item{col.point}{
  The color used for data points
  }
  \item{highlight}{
  The color used to higlight classes with mean sighnificantly different from zero
  }
  \item{standd}{
  Plot residuals against diameter standardized using the plot-specific mean and diameter (\code{standd=TRUE}) 
  or against raw diameter (\code{standd=FALSE})
  }
  \item{cex}{See \code{\link{par}}}
  \item{corD}{should predictions of random effects be plotted on mean diameter of the plot.}
  \item{ask}{ask before new plot.}
  \item{...}{Other arguments, currently ignored.}
}

\details{
The function makes residual plots on a fitted H-D model, which can be used to explore whether the fixed part 
satisfactorily models the shape of H-D models. 
The residuals are plotted on diameters standardized at plot level (dsd) or on raw diameters (d) according to argument \code{standd}. 
Here dsd=(d-D)/SDD, where d is tree diameter, D and SDD are the mean and standard deviation 
  of diameters on the plot in question.  Using plot-specific standardized diameter ensures that e.g.,
the medium-sized trees of the plot are always in the middle of the plot, which provides 
a better graph to explore the fit at the plot level in a dataset where the diameter range varies between plots.

}


\author{
Lauri Mehtatalo (lauri.mehtatalo@uef.fi)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The function plots model residuals on the required type of diameter and adds a whiskers plot 
using \code{\link[lmfor]{mywhiskers}} with argument \code{se=TRUE}.
}
\examples{
data(spati)

model<-fithd(spati$d,spati$h,spati$plot)

plot(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
