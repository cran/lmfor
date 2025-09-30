# TODO: Add comment
# 
# Author: 03191657
###############################################################################



# sama kuin ylla, mutta parametrit annetaan vektorina
# ja tuloksen attribuuttina palautetaan osittaisderivaatat parametrien suhteen
# thetan pitaa olla matriisi, jossa on n rivia ja sarakkeet kutakin primaarista 
# parametria varten
volvff<-function(dbh,h,theta=NA,logita=NA,lambda=NA,grad=TRUE) {
        if (!is.na(theta[1])) {
           logita<-theta[,1]
           lam<-exp(theta[,2])
           } else {
           lam<-exp(lambda)
           }
        w<-2-2*exp((h-1.3)/lam)/(1+exp((h-1.3)/lam))
        rstump<-w*dbh/20+(1-w)*dbh/20*h/(h-1.3)
        value<-pi*exp(logita)/(1+exp(logita))*(rstump)^2*(10*h)
        value<-pi*1/(1+exp(-logita))*(rstump)^2*(10*h)
        gr<-NA
        if (grad) { 
           gr<-cbind(logita=value/(1+exp(logita)),
                       lambda=-2*pi*dbh*h*rstump/(1+exp(-logita))*1.3/lam*
                       exp((h-1.3)/lam)/((1+exp((h-1.3)/lam))^2))
           }
        attr(value,"gradient")<-gr
        value
        }

# Laskee annetun datan puiden tilavuudet seka niihin 
# liittyvan parametrien estimoitivirheen vaikutuksen.
# palauttaa vektorin, jossa on 
# volsum - koknaistilavuus (motteina jos malli antaa ne litroina)
# vosumvar - koknaistilavuuden varianssi ja
# vosumci - kokonaistilavuuden 1-p % luottamusvali joka 
# Lisaksi attribuutteina palautetaan 
# trees - data.frame, jossa on jokaiselle puulle ennustettu tilavuus, 
#         sen odotusarvon estimaatin varianssiestimaatti ja 1-p% luottamsuvali 
# varmu - puukohtaisten tilavuusennusteiden varianssi-kovarianssimatriisi. 
# Luottamysvalin lasketaan 1. asteen taylor-approksimaation avulla linearisoidusta mallista.
# Ko. approksimaation avulla lasketaan ennustevektorin keskiarvon varianssi-kovarianssimatriisi ja sen alkiot summataan.  
# Kaikki varianssit  ja luottamusvalit siis ottavat huomioon vain 
# parametrien estimointivirheen vaikutuksen; eli jaannosvaihtelua ei ole huomioitu.   
predvff<-function(data,mod,p=0.05,
                  varMethod="taylor",
                  biasCorr="none",
                  nrep=500) {
         stopifnot(any(varMethod==c("taylor","simul","none")))
         stopifnot(any(biasCorr==c("none","integrate","twopoint")))
         data$logita<-data$lambda<-0
         n<-dim(data)[1]
         beta<-fixef(mod)
         varbeta<-mod$varFix
         formulas<-eval(mod$call$fixed)
         A1<-model.matrix(formulas[[1]],data=data)
         A2<-model.matrix(formulas[[2]],data=data)
         theta<-matrix(bdiag(A1,A2)%*%beta,nrow=n,byrow=FALSE)
         sdla<-mod$sigma*sqrt(sum(unlist(pdMatrix(mod$modelStruct$reStruct))))
         muhat0<-volvff(data$dbh,data$h,theta)
         if (biasCorr=="integrate") {
             muhat<-rep(NA,n)
             for (i in 1:n) {
             muhat[i]<-integrate(function(x) dnorm(x,mean=theta[i,1],sd=sdla)*
                                   volvff(data$dbh[i],data$h[i],cbind(x,theta[i,2])),
                                   qnorm(1e-10,mean=theta[i,1],sd=sdla),
                                   qnorm(1-1e-10,mean=theta[i,1],sd=sdla))$value
             }
             } else if (biasCorr=="twopoint") {
             muhat1<-volvff(data$dbh,data$h,theta+matrix(c(sdla,0),byrow=TRUE,nrow=n,ncol=2))
             muhat2<-volvff(data$dbh,data$h,theta-matrix(c(sdla,0),byrow=TRUE,nrow=n,ncol=2))
             muhat<-(muhat1+muhat2)/2
            } else { 
            muhat<-muhat0
            }
            
         if (varMethod=="taylor") {
            Xhat<-cbind(attr(muhat0,"gradient")[,1]*A1,attr(muhat0,"gradient")[,2]*A2)
            varmu<-Xhat%*%varbeta%*%t(Xhat)
            totvol<-sum(muhat)
            totvolvar<-sum(varmu)
            totvolci<-qnorm(c(p/2,1-p/2),mean=totvol,sd=sqrt(totvolvar))/1000
            qtree<-cbind(qnorm(p/2,mean=muhat,sd=sqrt(diag(varmu))),
                         qnorm(1-p/2,mean=muhat,sd=sqrt(diag(varmu))))
            } else if (varMethod=="simul")  {
            biascor<-muhat/muhat0
            betasim<-mvrnorm(nrep,mu=beta,Sigma=varbeta)
            thetasim<-bdiag(A1,A2)%*%t(betasim)
            musim<-matrix(ncol=nrep,nrow=n)
            cat("Simulating a total of",nrep,"replicates\n")
            cat("Now simulating replicate\n")
            for (i in 1:nrep) {
                cat("\b\b\b\b\b\b",i)
                thetathis<-matrix(thetasim[,i],nrow=n,byrow=FALSE)    
                musim[,i]<-biascor*volvff(data$dbh,data$h,thetathis)
                }
            varmu<-cov(t(musim))
            totvol<-sum(muhat)
            totvolsim<-apply(musim,2,sum)
            totvolci<-quantile(totvolsim,probs=c(p/2,1-p/2))/1e3
            totvolvar<-sum(varmu)
            qtree<-matrix(t(apply(musim,1,function(x) quantile(x,probs=c(p/2,1-p/2)))),byrow=FALSE,ncol=2)
            } else if (varMethod=="none") {
            totvol<-sum(muhat)
            totvolvar<-NA
            totvolci<-c(NA,NA)
            varmu<-matrix(NA,nrow=n,ncol=n)
            qtree<-matrix(NA,ncol=2,nrow=n)
            }
         value<-list(totvol=totvol/1000,totvolvar=totvolvar/1e6,totvolci=totvolci)
         attr(value,"trees")<-data.frame(volume=muhat,variance=diag(varmu),
                                         lb=qtree[,1],ub=qtree[,2])
         attr(value,"varmu")<-varmu
         value
         }

