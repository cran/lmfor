\name{NRnum}
\alias{NRnum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve a Systems of Nonlinear Equations Using the Newton's Method
}
\description{
Solves a system of equations of form \eqn{f_1(x) = 0},\eqn{f_2(x) = 0},...,\eqn{f_p(x) = 0} for 
vector \eqn{x} using the multidimensional version of the Newton-Raphson algorithm. 
The gradients  are solved numerically within the function using R-function \code{numericDeriv}. 
}
\usage{
NRnum(init, fnlist, crit = 6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
  vector of initial values for \eqn{x}.
  }
  \item{fnlist}{
  list of R-functions for \eqn{f_1(x)}, \eqn{f_2(x)}, ..., \eqn{f_p(x)} each function gets a vector-valued argument \eqn{x} and returns a scalar value.
  }
  \item{crit}{
  Convergence criterion. Stop iteration when (\eqn{|f_1(x)|+|f_2(x)|+...+|f_p(x)|<crit}). 
  }
  \item{\dots}{
  Other arguments passed to the functions of \code{fnlist}
  }

}

\value{
A list of components 
\item{par }{the value of vector \eqn{x} in the solution}
\item{crit }{the value of the convergence criterion at the solution}
If estimation fails (no solution is found during 100 iterations), both elements of the solution are NA's. 
}

\author{
Lauri Mehtatalo, <\email{lauri.mehtatalo@uef.fi}>
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

\keyword{ numerical computations }