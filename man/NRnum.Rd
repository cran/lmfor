\name{NRnum}
\alias{NRnum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve systems of equations using the Gauss-Newton algorithm
}
\description{
Solves systems of functions of form f_1(x) = 0, f_2(x) = 0,... for vector x using the Gauss-Newton algorithm
(the multidimensional version of the Newton-Raphson algorithm). The gradients 
are solved numerically within the function using R-function \code{numericDeriv}. 
}
\usage{
NRnum(init, fnlist, crit = 6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
  A vector of initial values for x.
  }
  \item{fnlist}{
  a list of R-functions for f_1(x), f_2(x), ... the functions get a vector-valued argument x and return a scalar value.
  }
  \item{crit}{
  The maximum accepted value of the convergence criteria. The applied criteria is the sum of absolute function values at the solution (|f_1(x)|+|f_2(x)+...|)
  }
  \item{\dots}{
  Other arguments passed to the functions of \code{fnlist}
  }

}

\value{
A list of components 
\item{par }{the value of vector x in the solution}
\item{crit }{the value of the convergence criterion at the solution}
If estimation fails (no solution is found during 100 iterations), both elements of the solution are NA's. 
}

\author{
Lauri Mehtatalo, lauri.mehtatalo@uef.fi
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Function \code{\link[lmfor]{NR}}.
}
\examples{
# Moment-based recovery of Weibull parameters
mu<-14
mu2<-210
muf<-function(theta) theta[2]*gamma(1+1/theta[1])-mu
mu2f<-function(theta) theta[2]^2*gamma(1+2/theta[1])-mu2
functions<-list(muf,mu2f)
momrec<-NRnum(c(3,13),functions)
momrec$par
}

