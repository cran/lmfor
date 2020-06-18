\name{fithd}
\alias{fithd}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a Height-Diameter model to forest tree data using functions of package \code{nlme}.
}

\description{
Fits either linear or nonlinear Height-Diameter (H-D) model into a dataset of tree heights and diameters. 
Possible hierarchy of the data can be taken into account through random effects. 
Several commonly used nonlinear two-parameter H-D functions are available. Linear functions can be used as well.  
}

\usage{
fithd(d, h, plot=c(), modelName="naslund", nranp=2,
      random=NA, varf=0, na.omit=TRUE, start=NA, bh=1.3,
	  control = list(), SubModels=NA, vfstart=0)
}

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
  
  \item{nranp, random}{
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

  \item{na.omit}{
  Should missing heights be omitted. Defaults to \code{TRUE}.
  }  
  
  \item{start}{
  A vector of the starting values of the parameters of the nlme fit. 
  If NA, then the starting values are computed using 
  the function computing the starting values (e.g., startHDnaslund, see \code{\link[lmfor]{HDmodels}}).}
  
  \item{bh}{
  The applied breast height. Defaults to 1.3 (meters).
  }
  
  \item{control}{
  Parameters to control of the model fitting algorithm, see \code{\link[nlme]{nlmeControl}} for details.
  }
  \item{SubModels}{Implemented only for nonlinear models. A character vector of length 2 or 3, 
  according to the number of parameters in the model. It allows submodels for parameters a, b (and c), where 
  the parameter is explaiend by plot-specific mean diameter ("~dmean"), plot-specific standard deviation "~dsd",
  or diameter standardized at plot level ("~dstd"), when the predictor is (d-D)/SDD (see teh documentation of argument varf). 
  Defaults to NA, which corresponds to no submodels, or \code{submodels=c("~1","~1","1")}}
  
  \item{vfstart}{
  Starting value of the power parameter delta of the variance function. Defaults to 0.
  }
}

\details{
Depending on the model (nonlinear or linear, mixed-effects model or marginal), the the model is 
fitted using one of the following functions functions of the \code{nlme} package:  
\code{nlme},  \code{lme}, \code{gls} or \code{gnls}. 

See available H-D functions at \code{\link[lmfor]{HDmodels}}. The user can define her own new functions
as specified at \code{\link[lmfor]{HDmodels}}.
}

\value{
An object of class \code{hdmod}, inheriting from class \code{nlme}. 
}
\references{
Mehtatalo, L., Gregoire, T.G., and de Miguel, S. Modeling Height-diameter curves for height prediction. Canadian Journal of Forest Research, 45(7): 826-837, \doi{10.1139/cjfr-2015-0054}
}
\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@uef.fi}>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[lmfor]{HDmodels}} for the available functions, Functions \code{\link{nlme}},  
\code{\link{lme}}, \code{\link{gls}} or \code{\link{gnls}} for details on model fitting, 
\code{\link[lmfor]{ImputeHeights}} for imputing unobserved tree heights.
}
\examples{
data(spati)

fithd(spati$d,spati$h,spati$plot)
fithd(spati$d,spati$h,spati$plot,SubModels=c("dmean","log(dmean)"),varf=2)
}

\keyword{ forest biometrics }
