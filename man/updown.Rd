\name{updown}
\alias{updown}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve a simple equation using a step halving algorithm.
}
\description{
Solves equations of form \eqn{f(x)=0}, for scalar \eqn{x} \eqn{(l<=x<=u)} using a simple step 
halving algorithm, where \eqn{f(x)} is a monotonic continuous function. 
Initial finite upper and lower bounds for 
x are required. The algorithm first computes \eqn{f} for \eqn{x=u} and \eqn{x=l}. 
If the sign was different then 
another call is performed at the midpoint \eqn{x=(u+l)/2}, and the midpoint is 
taken as a new upper or lower bound, according to the location of sign change. 
The upper or lower bound are repeatedly updated until the 
absolute value of f at the midpoint is below a specified criteria.    
}
\usage{
updown(l, u, fn, crit = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{
The initial lower bound
}
  \item{u}{
The initial upper bound
}
  \item{fn}{
R-function for \eqn{f(x)}
}
  \item{crit}{
The convergence criteria (Maximum accepted value of f at the solution is \eqn{10^{-\code{crit}}}).
}
}

\value{
A scalar giving the value of \eqn{x} at the solution. 
If the sign did not change between l and u, NA is returned. 
}

\author{
Lauri Mehtatalo <lauri.mehtatalo@uef.fi>
}
\section{Warning }{May lead to infinite loop for non-continuous functions. Works only with monotonic 
functions.} 

\examples{
## Compute the median of Weibull distibution
fn<-function(x) pweibull(x,5,15)-0.5
updown(1,50,fn)
}

\keyword{ numerical computations }