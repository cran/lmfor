\name{ImputeHeights}
\alias{ImputeHeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute missing tree heights into a forest data using a nonlinear (mixed-effects) model. 
}
\description{
A function to impute tree heights in a forest inventory situation where all trees have been measured 
for diameter but only some trees have been measured for height. 
}
\usage{
ImputeHeights(d, h, plot, modelName = "naslund", nranp = 2, varf = TRUE, 
              addResidual = FALSE, makeplot=TRUE, level = 1,
              start=NA, bh=1.3, control=list(),random=NA)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
  A numerical vector of tree diameters, usually given in cm.
  }

  \item{h}{
  A numerical vector of tree heights, usually given in meters. Should be of the same length as \code{d}. 
  }

  \item{plot}{
  A vctor of type \code{numeric} or \code{factor}, defining the groups of the data; usually the plot indices. 
  Should be of the same length as \code{d} and \code{h}.
  }

  \item{modelName}{
  Either (i) a character vector specifying the name of the nonlinear function or (ii) the formula specifying 
  a linear model.
  In case (i) the name should be one of the functions documented on the help page of \code{\link[lmfor]{HDmodels}}. 
  In case (ii), it should be the linear formula in the form that is entered to the function \code{\link[nlme]{lme}}, for example 
  \code{model=h~d+I(d^2)-1}
  }
  
  \item{nranp}{
  Parameters nranp and random specify two alternative ways to specify the random effects of the model. An easy but restricted way 
  is to use argument \code{nranp}. It is an integer between 0 and the number of fixed parameters, and has the the following meaning 
  in the case of nonlinear model: 
  \itemize{
     \item If \code{nranp=0}, then a model without random parameters is fitted. 
                              Results to a fixed-effects model, and argument \code{plot} is not used.  
     \item If \code{nranp=1}, then parameter a of a nonlinear function or the first coefficient of the linear formula 
                              is assumed to vary among plots. 
     \item If \code{nranp=2}, then a and b or the first two terms of the linear formula are assumed to vary among plots or  
     \item If \code{nranp=3}, then a b, and c of a three-parameter nonlinear model or three first coefficients of a linear model 
                              are assumed to vary among plots
     }
  In the case of linear model, the constant (if exists) it always counted as the first term.  
  
  As an alternative to nranp, argument \code{random} can be used to express the random part as a nlme formula, 
  but without specification of the grouping structure. The provided formula is passed to the lme or nlme 
  function. Argument random is always used when provided, so nranp has effect only if \code{random=NA} (the default).                           
  }
  

  \item{varf}{
  Numeric with values 0, 1 or 2. If 0 or FALSE, no variance function is used. 
  If varf=1, 2 or TRUE, then the power- type variance function var(e)=sigma^2*w^(2*delta) is used.
  where weight w is the raw diameter (when varf=1 or TRUE), or w=max(1,dsd+3) (when varf=2), 
  where dsd=(d-D)/SDD. Here d is tree diameter, D and SDD are the mean and standard deviation 
  of diameters on the plot in question.    
  }

  \item{addResidual}{
   Boolean. If \code{TRUE}, a random residual is added to the imputed height from a normal distribution using 
   the estimated variance function. If also \code{level=0} or if the plot did not include any measured heights
   to predict the random effects,  
   then also a randomly selected plot effect from among the predicted plot effects is added.
   The added plot effect is the same for all trees of a given plot.
   
}
  \item{makeplot}{
 Should a residual plot of the fitted model be produced for evaluation of goodness of fit? 
 The plot is produced using the default arguments of function \code{\link[lmfor]{plot.hdmod}}, 
 and is not affected by the value of arguments \code{level} and \code{addResidual}. 
}

  \item{level}{
The level of prediction. 0 means fixed-effect prediction and 1 means plot-level prediction using the random effects. Has no effect if \code{nranp=0}. 
}

  \item{start, bh, control, random}{
  Arguments passed to \code{fithd}. See documentation of \code{\link[lmfor]{fithd}}. 
}
}
\details{
The function predicts the missing 
heights using a nonlinear mixed-effects model or a nonlinear fixed-effects model. In mixed-effects model, 
plot-specific random effects can be used if other tree heights have been measured from the same plot. 
Also random, normally distributed residual can be added to the heights according to the estimated 
constant or heteroscedastic residual variance structure. 
}
\value{
A list of components
\item{h }{A vector of tree heights, including the measured heights for the trees with known height and imputed heights for the others.}
\item{imputed }{A booelan vector of the same length as h, having value TRUE for imputed heights. Produced as \code{is.na(data$h)}}
\item{model }{The fitted model that was used in imputation. Fitted using \code{\link[lmfor]{fithd}} which in turn 
calls \code{\link[nlme]{nlme}} or  \code{\link[nlme]{gnls}}}
\item{predType}{A vector of the same length as h, including information on the level of prediction. Value 0 means a measured height 
(no model prediction is used), value 1 means the plot-level prediction has been done using the estimated plot effects.
Value 2 means that no sample trees were available and the prediction is based on fixed part only (if level=0) or on a simulated 
plot effect (if level=1).}
\item{hpred}{Predicted heights for all trees. Equals to vector h for trees that had missing heights.}
}

\references{
Mehtatalo, L., de Miguel, S. and Gregoire, T.G. Modeling Height-diameter curves for prediction. Canadian Journal of Forest Research, 45(7): 826-837, 10.1139/cjfr-2015-0054.
}

\author{
Lauri Mehtatalo <lauri.mehtatalo@uef.fi>
}
\note{
Works only with the nonlinear functions specified in \code{\link[lmfor]{HDmodels}}; does not work if the modelName is specified as a linear expression.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[lmfor]{fithd}} for model fitting and \code{\link[lmfor]{plot.hdmod}} for plotting.
}
\examples{
data(spati) 

ImpFixed<-ImputeHeights(spati$d,spati$h,spati$plot,level=0)
ImpRandom<-ImputeHeights(spati$d,spati$h,spati$plot,level=1,makeplot=FALSE)
# Try also
# ImpRanRes<-ImputeHeights(spati$d,spati$h,spati$plot,level=1,addResidual=TRUE,makeplot=FALSE)

plot(spati$d[!is.na(spati$h)],
     spati$h[!is.na(spati$h)],
     col=spati$plot[!is.na(spati$h)],
     main="Observations", xlab="d, cm", ylab="h, m",
     ylim=c(0,30))

plot(spati$d[ImpFixed$imputed],
     ImpFixed$h[ImpFixed$imputed],
     col=spati$plot[ImpFixed$imputed],
     main="Imputed, Naslund, Fixed", xlab="d, cm", ylab="h, m",
     ylim=c(0,30))
       
plot(spati$d[ImpRandom$imputed],
     ImpRandom$h[ImpRandom$imputed],
     col=spati$plot[ImpRandom$imputed],
     main="Imputed, Naslund, Fixed + Plot", xlab="d, cm", ylab="h, m",
     ylim=c(0,30))

# Try also
# plot(spati$d[ImpRanRes$imputed],
#     ImpRanRes$h[ImpRanRes$imputed],
#     col=spati$plot[ImpRanRes$imputed],
#     main="Imputed, Naslund, Fixed + Plot + Tree", xlab="d, cm", ylab="h, m",
#     ylim=c(0,30))
}

\keyword{ forest biometrics}
