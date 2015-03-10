\name{recweib}
\alias{recweib}
\alias{func.recweib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recovery of Weibull parameters of tree diameter distribution using measured stand characteristics
}
\description{
The function finds such parameters shape and scale of the Weibull diameter distribution that 
yield the given basal area, number of stems and weighted/unweighted mean/median diameter. 
}
\usage{
recweib(G, N, D, Dtype, init=NA, trace=FALSE)
func.recweib(lshape, G, N, D, Dtype, trace=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{ The basal area in \eqn{m^2/ha}, scalar.  }
  \item{N}{The number of stems per ha, scalar. }
  \item{D }{Either A: The arithmetic mean diameter, B: The basal-area weighted mean diameter, 
      C: median diameter or D: The basal-area weighted median diameter of the stand, cm.}
\item{Dtype}{
One of characters "A", "B", "C", "D", indicating which type of mean diameter was given in argument D.
}
 \item{init}{The initial guess for the shape parameter (scalar). If not given, a simple model (see 
  Siipilehto and Mehtatalo 2013, appendix)
  is used to compute the initial guess. }

  \item{trace}{if TRUE, some output on the convergence of the algorithm is printed on the screen. 
}
\item{lshape}{logarithmic shape parameter, (log(shape+0.01))}
}

\details{
The recovery is based on the solution of the equation 
DQMW(shape,scale(D,shape))-DQM = 0, where DQMW(shape, scale(D,shape)) expresses the DQM of The Weibull distribution
for the given value of the shape parameter and using the scale parameter that corresponds 
to the given combination of the shape parameter and the mean/median diameter given in D. The function which is set to zero is implemented in \code{\link{func.recweib}}.
The Gauss-Newton method implemented in  \code{\link{NRnum}} is used 
for solving the equation. 
}
\value{
A list of components 
\item{shape, scale}{The value of the shape and scale parameters at the solution.}
\item{G, N, D, Dtype}{The input arguments.}
\item{val}{The value of the equation DQMW(shape,scale(D,shape))-DQM at the solution}
}
\references{
Siipilehto, J. and Mehtatalo, L. 2013. Parameter recovery vs. parameter prediction for the
Weibull distribution validated for Scots pine stands in
Finland. Silva Fennica 47(4), article id 1057.
}

\author{
Lauri Mehtatalo and Jouni Siipilehto
}

\seealso{
 The mean diameters for options A, B, C and D are computed by functions documented at 
\code{\link{scaleDMean}}.
}
\examples{
# Demonstration with 3 example stands.
# Example stand 1. Uneven-aged stand in Finland (Vesijako, Kailankulma, stand no 1):
G_ha<-17.0
N_ha<-1844
D<-7.9
DG<-19.6
DM<-8.1
DGM<-19.1
recweib(G_ha,N_ha,D,"A") #  1.066123,  8.099707
recweib(G_ha,N_ha,DG,"B") # 1.19316  8.799652
recweib(G_ha,N_ha,DM,"C") # 1.601795 10.18257
recweib(G_ha,N_ha,DGM,"D") # 1.095979 8.280063


# Example 2. Even aged stand in Finland (see Siipilehto & Mehtatalo, Fig 2):
G_ha<-9.6
N_ha<-949
D<-11.0
DG<-12.3
DM<-11.1
DGM<-12.4  
recweib(G_ha,N_ha,D,"A") # 4.465673 12.05919
recweib(G_ha,N_ha,DG,"B") # 4.463991 12.05912
recweib(G_ha,N_ha,DM,"C")  # 4.410773 12.05949
recweib(G_ha,N_ha,DGM,"D") # 4.448272 12.05924


# Example 3. Assumed peaked even aged stand (see Siipilehto & Mehtatalo, Fig 1):
G_ha<-10.0
N_ha<-1300
D<-9.89
DG<-10.0
DM<-9.89
DGM<-10.0  
recweib(G_ha,N_ha,D,"A")  #  34.542 10.04978
recweib(G_ha,N_ha,DG,"B") # 14.23261 10.22781
recweib(G_ha,N_ha,DM,"C") # 6.708882 10.44448
recweib(G_ha,N_ha,DGM,"D") # 24.45228 10.10607
}

\keyword{forest biometrics}
