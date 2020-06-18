\name{qqplotHD}
\alias{qqplotHD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal QQ-plot of a fitted H-D model
}
\description{
Produces a panel of graphs including the Normal qq-plot of a H-D model residuals and of the predicted random effects.
}
\usage{
qqplotHD(model, startnew=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A nonlinear H-D model model fitted by \code{\link[lmfor]{fithd}}.
  }
  
  \item{startnew}{
  Should a new plotting window be opened?
  }
}

\details{
The function extracts the residuals and the random effects of the fitted Height-Diameter model 
and produces a panel of plots including univariate Normal qq-plots of the model.
}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(spati)

model<-fithd(spati$d,spati$h,spati$plot)

qqplotHD(model)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ forest biometrics}
\keyword{ graphs} 
