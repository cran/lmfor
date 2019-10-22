\name{NR}
\alias{NR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve a Nonlinear Equation Using Newton-Raphson algorithm.
}
\description{
Solves an equation equation of form \eqn{f(x)=0}, for scalar \eqn{x} using the Newton-Raphson algorithm.  
}
\usage{
NR(init, fn, gr, crit = 6, range = c(-Inf, Inf))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
  Numeric scalar, The initial guess for \eqn{x}.
  }
  \item{fn}{
  An R-function returning the scalar value of \eqn{f(x)}, with \eqn{x} as the only argument. 
  }
  \item{gr}{
  An R-function returning the first derivative of \eqn{f(x)}, with \eqn{x} as the only argument. 
  }
  \item{crit}{
  Convergence criteria. The upper limit for the absolute value of \eqn{f(x)} at an accepted the solution. 
  }
  \item{range}{
  A two-unit vector giving the upper and lower bounds for \eqn{x}. The solution is searched from within this range.   
  }
}

\details{
The function is a straightforward implementation of the well-known Newton-Raphson algorithm. 
}

\value{
A list of components 
\item{par }{the value of \eqn{x} in the solution}
\item{crit }{the value of \eqn{f(x)} at the solution}
If estimation fails (no solution is found during 100000 iterations), both e
lements of the solution are NA's. 
}
\author{
Lauri Mehtatalo <lauri.mehtatalo@uef.fi>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link[lmfor]{NRnum}} for a vector-valued x without analytical gradients.
}
\examples{
## Numerically solve Weibull shape for a stand
## where DGM=15cm, G=15m^2/ha and N=1000 trees per ha
func<-function(shape,G,N,DGM) {
##      print(G,DGM,N)
      val<-pi/(4*gamma(1-2/shape)*log(2)^(2/shape))-G/(N*DGM^2)
      val
      }

grad<-function(shape) {
      pi/4*(-1)*                                                            
      (gamma(1-2/shape)*log(2)^(2/shape))^(-2)*                             
      (gamma(1-2/shape)*digamma(1-2/shape)*2*shape^(-2)*log(2)^(2/shape)+   
      log(2)^(2/shape)*log(log(2))*(-2)*shape^(-2)*gamma(1-2/shape))         
      }
      
shape<-NR(5,fn=function(x) func(x,G=10000*15,1000,15),gr=grad,crit=10,range=c(2.1,Inf))$par
}

\keyword{ numerical computations }
% Add one or more standard keywords, see file 'KEYWORDS' in the

