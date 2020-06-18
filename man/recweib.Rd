\name{recweib}
\alias{recweib}
\alias{func.recweib1}
\alias{func.recweib2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recovery of Weibull parameters of tree diameter distribution using measured stand characteristics
}
\description{
The function finds such parameters shape and scale of the Weibull diameter distribution that 
yield the given basal area, number of stems and weighted/unweighted mean/median diameter. 
Weibull function can be assumed either as the unweighted or basal-area weighted distribution. 
}
\usage{
recweib(G, N, D, Dtype, init=NA, trace=FALSE, weight=0)
func.recweib1(lshape, G, N, D, Dtype, trace=FALSE)
func.recweib2(lshape, G, N, D, Dtype, trace=FALSE)
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
  is used to compute the initial guess for the unweighted case; value 4 is used as default in the basal-area weighted case. }

  \item{trace}{if TRUE, some output on the convergence of the algorithm is printed on the screen. 
}
  \item{weight}{if weight=0 (the default), Weibull function is assumed as the unweighted density. If weight=2, weibull 
                function is assumed as the basal.area weighted density. 
}
\item{lshape}{logarithmic shape parameter, (log(shape+0.01))}
}

\details{
The recovery is based on the solution of the equation 
DQMW^2(shape,scale(D,shape))-DQM^2= 0, where DQMW(shape, scale(D,shape)) expresses the DQM of the assumed Weibull distribution
for the given value of the shape parameter and using the scale parameter that corresponds 
to the given combination of the shape parameter and the mean/median diameter given in D. 
The function which is set to zero is implemented in functions \code{\link{func.recweib1}} (unweighted case) and 
\code{\link{func.recweib2}} (ba-weighted case).
The Gauss-Newton method implemented in  \code{\link{NRnum}} is used 
for solving the equation. 
}
\value{
A list of components 
\item{shape, scale}{The value of the shape and scale parameters at the solution.}
\item{G, N, D, Dtype}{The input arguments.}
\item{val}{The value of the equation DQMW^2(shape,scale(D,shape))-DQM^2 at the solution}
}
\references{
Siipilehto, J. and Mehtatalo, L. 2013. Parameter recovery vs. parameter prediction for the
Weibull distribution validated for Scots pine stands in
Finland. Silva Fennica 47(4), article id 1057. \doi{10.14214/sf.1057}

Mehtatalo, Lauri and Lappi, Juha 2020a. Biometry for Forestry and Environmental Data: with examples in R. New York: Chapman and Hall/CRC. 426 p. \doi{10.1201/9780429173462}

Mehtatalo, Lauri and Lappi, Juha 2020b. Biometry for Forestry and Environmental Data: with examples in R. Full Versions of The Web Examples. Available at \url{http://www.biombook.org}. 

}

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}> and Jouni Siipilehto
}

\seealso{
 The mean diameters for options A, B, C and D are computed by functions documented at 
\code{\link{scaleDMean1}}.
}
\examples{
# Demonstration with 3 example stands.
# Example stand 1. Uneven-aged stand in Finland (Vesijako, Kailankulma, stand no 1):
G<-17.0
N<-1844
D<-7.9
DG<-19.6
DM<-8.1
DGM<-19.1
recweib(G,N,D,"A")            #  1.066123 8.099707
recweib(G,N,DG,"B")           # 1.19316  8.799652
recweib(G,N,DM,"C")           # 1.601795 10.18257
recweib(G,N,DGM,"D")          # 1.095979 8.280063
recweib(G,N,D,"A",weight=2)   # 2.590354 21.66398
recweib(G,N,DG,"B",weight=2)  # 2.563892 22.07575
recweib(G,N,DM,"C",weight=2)  # 2.998385 17.74291
recweib(G,N,DGM,"D",weight=2) # 2.566621 22.03183



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
\keyword{diameter distributions}
