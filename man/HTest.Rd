\name{HTest}
\alias{HTest}
\alias{area_esh}
\alias{gg_wind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate stand density using a Horvitz--Thompson-like estimator
}
\description{
\code{HTest} calculates the Horvitz--Thompson-like stand density estimate (number of trees) in a specified area based on a collection of detected trees.

\code{area_esh} is an internal function for surface area calculations that can handle empty sets.

\code{gg_wind} is an internal function that forms a union of discs based on their center points and radii.
}
\usage{
HTest(treelist, plotwindow, alpha)

area_esh(W)

gg_wind(treelist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treelist}{
A 3-column matrix containing the x and y coordinates of detected trees and their crown radii.
}
   \item{plotwindow}{
A \code{spatstat} object of class ''\code{owin}'', representing the area where stand density estimation is done.
}
 \item{alpha}{
A tuning parameter that controls the calculation of detection probabilities, or detectabilities. Must have a value from -1 to 1.
}

 \item{W}{
A \code{spatstat} object of class ''\code{owin}'' or NULL.
}
}
\details{
\code{HTest} is the Horvitz--Thompson-like stand density estimator presented by Kansanen et al. (2016) to adjust individually detected trees for non-detection. It uses individual tree detection data, namely the locations and crown radii of detected trees, to calculate detection probabilities, or detectabilities, for every detected tree, and produces an estimate based on the detectabilities. The detectability for a certain tree is based on the planar set formed by the larger trees. The parameter \code{alpha} controls how easy it is to detect a tree of certain size from under the larger trees. If \code{alpha}=1, then the tree will be detected if it is not fully covered by the larger crowns. If \code{alpha}=0, the tree will be detected if its center point is not covered. If \code{alpha}=-1, the tree will be detected if it is fully outside the larger tree crowns.

The object \code{treelist} can include trees that are not in the estimation area specified by \code{plotwindow}. This can be useful to take into account possible edge effects, by including trees with center points outside \code{plotwindow} that have crown discs that intersect \code{plotwindow}. The estimate is calculated only using those trees that have crown center points in \code{plotwindow}.

\code{area_esh} and \code{gg_wind} are internal helper functions used by \code{HTest}. First one is a shell for the \code{spatstat} function \code{area.owin}  that takes into account that an intersection of two sets can be empty, represented in the calculations as NULL. The function returns 0 in this case. Otherwise, it returns the surface area of the window \code{W}. The latter function forms a union of discs that is needed in the detectability calculations. 

}
\value{
\code{HTest} returns a list with two components:
\item{N }{The estimated number of trees in \code{plotwindow}}
\item{treelist }{matrix with columns ''r'' and ''detectability'', giving the tree crown radii that have been used in the estimation, as well as the detectabilities for trees with those crown radii.}

\code{area_esh} returns 0, if \code{W} is NULL; otherwise, the surface area of \code{W}.

\code{gg_wind} returns a \code{spatstat} object of class ''\code{owin}'' representing a set formed as a union of discs.
}
\references{
Kansanen, K., Vauhkonen, J.,  Lahivaara, T., and Mehtatalo., L. (2016) \emph{Stand density estimators
based on individual tree detection and stochastic geometry}. Canadian Journal of Forest Research 46(11):1359--1366. http://dx.doi.org/10.1139/cjfr-2016-0181.

Baddeley, A., Rubak, E. and Turner, R. (2015) \emph{Spatial Point Patterns: Methodology and Applications with {R}}. Chapman and Hall/CRC Press, London. http://www.crcpress.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/9781482210200/

Kansanen, K., Packalen, P., Lahivaara, T., Seppanen, A., Vauhkonen, J., Maltamo, M.,  and Mehtatalo., L. (2019) \emph{Horvitz--Thompson-like stand density estimation and functional k-NN in individual tree detection}. Submitted manuscript.
}
\author{
Kasper Kansanen, kasperkansanen@gmail.com
}
\note{
These functions require the package \code{spatstat} (Baddeley et al. 2015) to work.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
require(spatstat)

# Generate a 10x10 meter square window:
w<-square(10)

# Generate 6 detected trees, 5 located in the window:
x<-cbind(c(6.75, 8.65, 3.95, 2, 2, 11), 
         c(1.36, 3.10, 6.66, 2, 4, 11), 
         c(1.29, 2.31, 1.80, 2, 1.5, 3))

# Draw the set formed by  the detected tree crowns:
plot(w)
plot(gg_wind(x), add=TRUE)

# Calculate the results with different alpha:
HTest(x, w, 1)
HTest(x, w, 0)
HTest(x, w, -0.75)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
