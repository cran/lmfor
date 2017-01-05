\name{mywhiskers}
\alias{mywhiskers}
\title{A whiskers type residual plot}
\description{
A function for adding vertical lines onto residual plots to show 
\enumerate{
\item 95\% confidence intervals of means or 
\item 95\% confidence intervals for individual observations 
}in the classes of the 
variable on the x-axis. Plot of the first type is useful for analyzing the fit of the assumed fixed part and plots of 
type b can be used to analyze the homogeneity of residuals. 
}
\usage{
mywhiskers(x, y, 
           nclass = 10, 
           limits = NA, 
           add = FALSE, 
           se = TRUE, 
           main = "", 
           xlab = "x", 
           ylab = "y", 
           ylim = NA, 
           lwd = 1, 
           highlight = "red")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The variable on the x-axis. Usually one of the predictors or the predicted value.
}
  \item{y}{The variable on the y-axis. Usually model residual.
}
  \item{nclass}{
The maximum number of classes to be used.
}
  \item{limits}{
The class limits. Alternative to nclass.
}
  \item{add}{
logical. Whether a new graphic window is opened or the lines will be added into an exosting plot.
}
  \item{se}{
Logical. Use standard errors of means (se=TRUE, option (a) above) or class-specific standard deviations 
(se=FALSE, option (b) above).
}
  \item{main, xlab, ylab, ylim, lwd}{
Graphical parameters of the plot. ignored if \code{add=TRUE}.
}
  \item{highlight}{
The color for lines that do not cross the y-axis.
}
}
\details{
The function first classifies the data in \code{nclass} classes of variable x so that 
each class has approximately equal number of observations. 
Then the class mean and deviation s is computed for each class, where s is either 
the standard error of the mean (if \code{se=TRUE}) 
or standard deviation (if \code{se=FALSE}). 
A vertical line is plotted at the middle of each class showing the class mean by a dot and 
lines of length 3.92*s. If the line does not cross the x- axis, then the highlight color is used in the line.
With small number of observations (or lot of ties), the number of classes is decreased until 
each class includes the minimum of 2 observations.  
}
\value{
The function is usually used for its side effects (i.e., the plot). However, the values used in producing 
the plot are returned in a list of elements 
\itemize{
\item x: the class middlepoint \code{x} values.
\item m: class-specific means of \code{y}.
\item s: class-specific standard deviations or standard errors of \code{y} (see details).
\item lb: lower ends of the class-specific lines.
\item ub: upper ends of the lines.
}
}
\references{
Mehtatalo, L. 2010. Forest biometrics with R. Lecture notes.
}
\author{
Lauri Mehtatalo
}

\examples{
x<-seq(1,100,1)
y<-x+10*log(x)+rnorm(100,0,5)
fm1<-lm(y~x)
plot(x,resid(fm1))
mywhiskers(x,resid(fm1),se=FALSE,add=TRUE)
mywhiskers(x,resid(fm1),se=TRUE,lwd=2,add=TRUE)
abline(h=0)
}

\keyword{ graphs} 