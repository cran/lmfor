\name{HTest_cps}
\alias{HTest_cps}
\alias{detectability_cps}
\alias{visibility_thinning_cps}
\alias{ordering_cps}
\alias{polar_to_cart}
\alias{cart_to_polar}
\alias{triangle_coords}
\alias{shades}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate forest characteristics of interest in circular plot sampling using a Horvitz--Thompson-like estimator
}
\description{
\code{HTest_cps} calculates Horvitz--Thompson-like estimates of forest characteristics of interest in a specified circular area based on a collection of detected trees and their detection probabilities, or detectabilities. Also produces estimated variances and confidence intervals.

\code{detectability_cps} calculates detectabilities of trees in a circular plot sample.

\code{visibility_thinning_cps} takes a tree list and determines if the trees can be detected when a certain visibility-based detection condition is used.

\code{ordering_cps} is a helper function for preprocessing of tree lists: it takes a tree list and orders the trees based on their distance to plot centre point.

\code{polar_to_cart} and \code{cart_to_polar} are internal functions for transforming polar coordinates to cartesian coordinates and vice versa.

\code{triangle_coords} is an internal function that, given locations and diameters of discs, returns coordinates needed to define the areas behind the discs that are non-visible from the origin.

\code{shades} is an internal function that forms polygonal approximations of the planar sets that are non-visible from the origin.
}
\usage{
HTest_cps(data, total=TRUE, confidence.level=0.95)
detectability_cps(data, plot.radius, alpha=0, polar=TRUE, npoly=1024, delta=NULL)
visibility_thinning_cps(data, plot.radius, alpha=0, polar=TRUE, npoly=1024, delta=NULL)
ordering_cps(data, polar=TRUE)
polar_to_cart(X)
cart_to_polar(X)
triangle_coords(X, plot.radius, polar=TRUE)
shades(X, plot.radius, polar=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
For \code{HTest_cps} a matrix where each row corresponds to a tree, the first column contains the detectabilities of the trees calculated with \code{detectability_cps}, and other columns correspond to measurements over which estimation is wanted. For \code{detectability_cps} a four-column matrix, each row containing the coordinates of a tree, the diameter of the tree, and an indicator if the tree has been detected (1) or not (0). In other words, each row is a vector of the form (r, phi, diameter, detected) if the locations are given in polar coordinates, or (x, y, diameter, detected) if the locations are given in cartesian coordinates. It is assumed that the trees are a sample from circular plot sampling, and that the centre point of the plot is the origin (0,0). For \code{visibility_thinning_cps} a three-column matrix, otherwise similar to the matrix needed by \code{detectability_cps} but without the ''detected'' column.
For \code{ordering_cps} a matrix, each row corresponding to a tree, the first three columns being either (x, y, diameter) or  (r, phi, diameter).
}

\item{total}{
Do you want to estimate population totals (TRUE) or population means (FALSE)?
}

\item{confidence.level}{
The level of the approximate confidence interval, a value between 0 and 1. For example, \code{confidence.level}=0.95 indicates a 95 per cent confidence interval. 
}

\item{plot.radius}{
Radius of the plot in which circular plot sampling has been performed.
}

\item{alpha}{
A tuning parameter that controls the calculation of detection probabilities, or detectabilities.
}

\item{polar}{
Are the locations of trees given in polar coordinates (TRUE) or cartesian coordinates (FALSE)?
}

\item{npoly}{
Number of edges for the polygonal approximation of the plot boundary and the circles used to calculate the detection probabilities. Used if \code{delta}is NULL.
}

\item{delta}{
The tolerance of the polygonal approximation of of the plot boundary and the circles used to calculate the detection probabilities: the length of the arc that will be replaced by one edge of the polygon. If given value that is different from NULL tihs will override \code{npoly}
}

\item{X}{
For \code{polar_to_cart} a two-column matrix, each row containing the polar coordinates (r, phi) of a point. For \code{cart_to_polar} a two-column matrix, each row containing the cartesian coordinates (x, y) of a point. For \code{triangle_coords} and \code{shades} a three-column matrix, each row containing either the polar coordinates and a diameter of a tree (r, phi, diameter) or the cartesian coordinates and a diameter of a tree (x, y, diameter). The coordinates define the centre point of the tree.
}

}
\details{
The function \code{HTest_cps} produces estimates of forest characteristics of interest in a circular plot sampling situation. More specifically, it is assumed that an observer stands in a point and observes such trees that are within the fixed-area plot and are not hiding behind other trees. It is assumed that observer can record the locations and diameters of the trees that they observe. The observer can be a person or a piece of equipment, such as terrestrial laser scanner or camera.

The estimation is based on a Horvitz--Thompson-like estimator presented by Kansanen et al. (2019). This construction uses approximated detection probabilities, or detectabilities, that depend on the size and distance from the plot centre of the tree for which the probability is calculated, the nonvisible area produced by trees that are closer to the centre point, and a visibility based detection condition.  It is assumed that the centre point of the sampling plot is the origin of the plane (0,0). The function \code{detectability_cps} is used to calculate the detectabilities.

The confidence intervals that \code{HTest_cps} produces are based on the t-distribution if less than 50 trees have been observed, and the standard normal distribution otherwise.

The parameter \code{alpha} is a value between -1 and 1 and it controls the detection condition. alpha=1 means that trees are detected if the stems are fully visible to the observer, alpha=0 means that they are detected if the center point is visible, and alpha=-1 means that a tree is detected if any part of the stem is visible.

The estimation is not possible if the data contains trees that cover the plot centre point. 

All of the variables related to distance and size, meaning the cartesian coordinates, the distance coordinate, \code{plot.radius}, \code{delta}, and tree diameters, should have the same unit, e.g. they should all be in metres.

The function \code{visibility_thinning_cps} is useful for simulation-based testing of the estimator. Given a tree list, it classifies trees as either detected or not detected based on a visibility based detection condition.

The function \code{ordering_cps} is a useful preprocessing step for tree lists over which estimation is needed. It reorders the rows of the tree list, corresponding to trees, based on the closest distance from the stem disc to the plot centre point. \code{detectability_cps} and \code{visibility_thinning_cps} use this ordering for their advantage, as this is the assumed order or sequence of detection. Be warned that even if you do not order the data with \code{ordering_cps} these functions will, and will output tree lists with this ordering!

The functions \code{polar_to_cart}, \code{cart_to_polar}, \code{triangle_coords}, and \code{shades} are internal functions used by the two main functions. They can be useful for visualizing data.
}
\value{
\code{HTest_cps} returns a four-column matrix, the columns containing an estimate, estimated variance of the estimator, and lower and upper bound of the approximate confidence interval for the estimate. Rows of the matrix correspond to the forest characteristics of interest given in columns \code{2:ncol(data)} of the input matrix \code{data}. If the input matrix has named columns, these names are used as row names of the output matrix.

\code{detectability_cps} returns a matrix with the locations and diameters of the trees given as input, the indicators of their detection, and the estimated detection probabilities.

\code{visibility_thinning_cps} returns a four-column matrix with the locations and diameters of the trees given as input and the indicators of their detection: 1, if a tree has been detected based on the detection condition given by \code{alpha}, and 0, otherwise.

\code{ordering_cps} returns a matrix with same dimensions as the input matrix, rows being reordered in the assumed order of detection.

\code{polar_to_cart} returns a two-column matrix of cartesian coordinates (x, y).

\code{cart_to_polar} returns a two-column matrix of polar coordinates (r, phi).

\code{triangle_coords} returns an eight-column matrix, each row containing cartesian coordinates needed for forming a polygonal representation of an area behind a tree, nonvisible from the origin.

\code{shades} returns a list of \code{owin} objects, each representing an area behind a tree, nonvisible from the origin.
}
\references{
Baddeley, A., Rubak, E. and Turner, R. (2015) \emph{Spatial Point Patterns: Methodology and Applications with {R}}. Chapman and Hall/CRC Press, London. http://www.crcpress.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/9781482210200/

Kansanen, K., Packalen, P., Maltamo, M., and Mehtatalo, L. (2019) \emph{Horvitz--Thompson-like estimation with distance-based detection probabilities for circular plot sampling of forests}. Submitted manuscript.

}
\author{
Kasper Kansanen, kasperkansanen@gmail.com
}
\note{
These functions require the package \code{spatstat} (Baddeley et al. 2015) to work.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Simulate some data from the Poisson process, in a plot of radius 10 metres 
# and stem density of 1000 trees/ha:
N <- rpois(1, lambda=0.05*pi*10^2)
proc <- cbind(10*sqrt(runif(N, 0 ,1)),  runif(N, -pi, pi), rweibull(N, shape=12, scale=22)/100)

# Preprocess the data to the right order:
proc <- ordering_cps(proc)

# Thin the data:
thinned<-visibility_thinning_cps(data=proc, plot.radius=10, alpha=1)

# Calculate the detection probabilities:
detprob <- detectability_cps(data=thinned, plot.radius=10, alpha=1)

# Calculate estimates of stem density (number of trees) and basal area 
# (the sum of areas covered by tree stems at breast height):
mydata<-cbind(detprob[,5], 1, pi*detprob[,3]^2)

HTest_cps(data=mydata)

# Calculate estimate of mean diameter at breast height with a 99 per cent 
# approximate confidence interval:
mydata<-cbind(detprob[,5], detprob[,3])

HTest_cps(data=mydata, total=FALSE, confidence.level=0.99)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

