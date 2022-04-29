\name{volvff}
\alias{volvff}
\alias{predvff}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Variable Form-Factor Volume Model
}
\description{
An R-function for the variable form-factor volume model and a function for computing 
bias-corrected volumes augmented with parameter uncertainty from a model fitted using nlme.}
\usage{volvff(dbh,h,theta=NA,logita=NA,lambda=NA)
predvff(data,mod,p=0.05,varMethod="taylor",biasCorr="none",nrep=500)
}

\arguments{
  \item{dbh}{
  vector of individual tree diameters, cm
  }
  \item{h}{
  vector of individual tree heights (m), of same length as \code{dbh}.
  }
  \item{theta, logita, lambda}{
  The parameters \code{logita} and \code{lambda} can be given either in a two-column 
  matrix where the number of rows equals to the length of \code{dbh}, 
  or separately in objects \code{logita} and \code{lambda}, which can be 
  scalars or vectors of same length as \code{dbh}. 
  If \code{theta} is not \code{NA}, argumets \code{logita} 
  and \code{lambda} are ignored.  
  }
  \item{data}{
  A data set including variables \code{dbh} and \code{h} and 
  all variables that are used as second-order predictors in \code{mod}. 
  }
  \item{mod}{
  A variable form-factor model fitted using \code{nlme}. 
  Any second-order predictors can be used for parameters \code{logita} and 
  \code{lambda}. The bias-correction works properly only 
  if the model includes random intercepts in \code{logita}.   
  }
  \item{p}{
  The probability used in constructiong the confidence intervals for 
  tree-level volumes and total volume. 
  Symmetric \eqn{100(1-p)}\% confidence intervals are returned.  
  }
  \item{varMethod}{
  Either \code{"taylor"}(the default) or \code{"simul"}. Specifies whether parameter 
  uncertainty is approximated by using a linearized model based on the 
  first-order Taylor approximation or by using Monte Carlo simulation.  
  }
  \item{biasCorr}{
  Either \code{"none"}(the default), \code{"integrate"} or \code{"twopoint"}. 
  Specifies whether the bias correction is not done (\code{"none"}), 
  is done by computing expected value of the prediction over the distribution 
  of random intercepts in parameter \code{logita}, or by using the two-point 
  approximation described in Kangas et al. (XXXX).  
  }
  \item{nrep}{
  The number of replicates in the Monte Carlo simulation when 
  \code{varMethod="simul"}. Ignored if \code{varMethod="taylor"}.
  }
  }

\details{
The variance-form-factor function is of form
\deqn{
v(D,H,a,\lambda)=\pi \frac{\exp(a)}{1+\exp(a)} R(D,H,\lambda)^2 H
}
where \eqn{a} is the logit-transformed form factor and \eqn{R(D,H,\lambda)} 
is the stem radius at stump height, which is approximated using 
\deqn{
R(D,H,\lambda)=w(H,\lambda)\frac{D}{2}+(1-w(H,\lambda))\frac{H}{H-B}\frac{D}{2}
}
where the weight is taken from the right tail of the logit transformation 
\deqn{
    w(H,\lambda)=2-2\frac{\exp\left( \frac{H-B}{\exp(\lambda)}\right)}{1+\exp\left( \frac{H-B}{\exp(\lambda)}\right)}
    }

Parameter uncertainty is reported because the same realized errors are used always when 
a model based on certain model fitting data is used; therefore those errors behave 
in practice like bias. Variance for total volume is computed as sum of all elements of the 
variance-covariance matrix of prediction errors of the mean. 
}

\value{
Function \code{volvff} returns a vector of tree volumes (in liters) that is of same length as 
vector \code{dbh}. In addition, attribute \code{grad} returns the Jacobian, which is used 
in nlme fitting for computing the derivatives of the model with respect to parameters, 
and in approximating the parameter uncertainty when \code{varMethod="taylor"}. 

Function \code{predvff} returns a list with following objects
  \item{totvol}{
  Total volume of the trees of \code{dat}, in \eqn{m^3}.
  }
  \item{totvolvar}{
  The estimated variance of \code{totvol}, taking into account parameter uncertainty. 
  }
  \item{totvolci}{
  The estimated \eqn{100(1-p)}\% confidence interval of \code{totvol}, based on parameter uncertainty. 
  If \code{varMethod="taylor"}, it is based on normality and 
  the approximated variance \code{totvolvar}. 
  If \code{varMethod="simul"}, the empirical quantiles of the 
  \code{nrep} simulations are used. 
  }
And attributes
  \item{trees}{
  A data frame of including tree-level volumes, their variance and 95\% confidence intervals.  
  }
  \item{varmu}{
  The variance-covariance matrix of prediction errors, taking into account theparameter uncertainty. 
  }
  }

\author{
Lauri Mehtatalo <\email{lauri.mehtatalo@luke.fi}>
}
\examples{
\dontrun{
library(lmfor)
data(treevol)

treevol$formfactor<-treevol$v/volvff(treevol$dbh,treevol$h,logita=100,lambda=log(0.2))
treevol$logitff<-log((treevol$formfactor)/(1-(treevol$formfactor)))

ptrees<-treevol[treevol$species=="pine",]
mod.init<-lm(logitff~I(1/h)+h+dbh+I(h*dbh)+I(1/(h*dbh))+
                     dataset+dataset:dbh+dataset:h,data=ptrees)

mod<-nlme(v~volvff(dbh,h,logita=logita,lambda=lambda),
         fixed=list(logita~I(1/h)+h+dbh+I(h*dbh)+I(1/(h*dbh))+dataset+
                           dataset:dbh+dataset:h+soil+temp_sum,
                           lambda~1),
         random=logita~1|stand/plot,
         start=c(coef(mod.init),rep(0,2),log(0.2)),
         data=ptrees,
         weights=varComb(varIdent(form=~1 |dataset),varPower()),
         method="ML",
#         control=list(msVerbose=TRUE),
#         verbose=TRUE
)

pred1<-predvff(ptrees,mod,varMethod="simul",biasCorr="integrate") 
pred1$totvol    
pred1$totvolvar
pred1$totvolci         
head(attributes(pred1)$trees)    

pred2<-predvff(ptrees,mod,varMethod="taylor",biasCorr="twopoint") 
pred2$totvol    
pred2$totvolvar
pred2$totvolci         
head(attributes(pred2)$trees)   
}
}
\references{
Kangas A., Pitkanen T., Mehtatalo L., Heikkinen J. (xxxx)
Mixed linear and non-linear tree volume models with regionally varying parameters
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree volume}
\keyword{forest biometrics}
