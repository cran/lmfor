# TODO: Add comment
# 
# Author: lamehtat
###############################################################################

dll<-function(x,mu,sigma,xi=0,lambda=1,log=FALSE) {
#	sigma<-exp(lsigma)/(1+exp(lsigma))
#	xi<-exp(lxi)
#	lambda<-exp(llambda)
	f<-lambda/sigma*
			1/((x-xi)*(xi+lambda-x))*
			1/(exp(-mu/sigma)*((x-xi)/(xi+lambda-x))^(1/sigma)+
				exp(mu/sigma)*((x-xi)/(xi+lambda-x))^(-1/sigma)+2)
	f[x<=xi|x>=(xi+lambda)]<-0
	f[lambda<=0|sigma<=0|sigma>=1]<-NaN
	if (log) f<-log(f)
	f
}

# cdf of logit-logistic distribution     
pll<-function(q,mu,sigma,xi=0,lambda=1,lower.tail=TRUE,log.p=FALSE) {
#	sigma<-exp(sigma)/(1+exp(sigma))
	F<-1/(1+exp(mu/sigma)*((q-xi)/(xi+lambda-q))^(-1/sigma))
	F[q<=xi]<-0
	F[q>=xi+lambda]<-1
	F[lambda<=0|sigma<=0|sigma>=1]<-NaN
	if (!lower.tail) F<-1-F
	if (log.p) F<-log(F)
	F
}

qll<-function(p,mu,sigma,xi=0,lambda=1,lower.tail=TRUE,log.p=FALSE) {
    if (log.p) p<-exp(p)
	if (!lower.tail) p<-1-p
#	sigma<-exp(sigma)/(1+exp(sigma))
	a<-((1/p-1)/exp(mu/sigma))^(-sigma)
	value<-(xi+a*xi+a*lambda)/(1+a)
	value[p==1]<-xi+lambda
	value[p<0|p>1]<-NaN
	value[lambda<=0|sigma<=0|sigma>=1]<-NaN
	value
}

# random number generation from the logit-logistic distribution         
rll<-function(n,mu,sigma,xi=0,lambda=1) {
     if (length(n)>1) n<-length(n)
     xi<-rep(xi,length.out=n)
	 lambda<-rep(lambda,length.out=n)
	 sigma<-rep(sigma,length.out=n)
	 mu<-rep(mu,length.out=n)
     qll(runif(n),mu,sigma,xi,lambda)
     }


pPercbas<-function(q,xi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
         if (length(xi)!=length(F)) stop("F and xi should be of equal length")
         if (length(xi)<2) stop("The number of percentiles should be >1")
         k<-length(F)
         if (min(xi[-1]-xi[-k])<=0) stop("The elements in xi should be strictly increasing")
         if (min(F[-1]-F[-k])<=0) stop("The elements in F should be strictly increasing")
         if (F[1]!=0|F[k]!=1) stop("F[1] should be 0 and F[k] should be 1")
         approx(xi,F,xout=q,yleft=0,yright=1)$y
         }

qPercbas<-function(p,xi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
          if (length(xi)!=length(F)) stop("F and xi should be of equal length")
          if (length(xi)<2) stop("The number of percentiles should be >1")
          k<-length(F)
          if (min(xi[-1]-xi[-k])<=0) stop("The elements in xi should be strictly increasing")
          if (min(F[-1]-F[-k])<=0) stop("The elements in F should be strictly increasing")
          if (F[1]!=0|F[k]!=1) stop("F[1] should be 0 and F[k] should be 1")
		  if (min(p)<0|max(p)>1) stop("p should be within [0,1]")
          approx(F,xi,xout=p,yleft=NaN,yright=NaN)$y
}

dPercbas<-function(x,xi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
	if (length(xi)!=length(F)) stop("F and xi should be of equal length")
	if (length(xi)<2) stop("The number of percentiles should be >1")
	k<-length(F)
	if (min(xi[-1]-xi[-k])<=0) stop("The elements in xi should be strictly increasing")
	if (min(F[-1]-F[-k])<=0) stop("The elements in F should be strictly increasing")
	if (F[1]!=0|F[k]!=1) stop("F[1] should be 0 and F[k] should be 1")
	b<-c(0,(F[-1]-F[-k])/(xi[-1]-xi[-k]),0)
	xi2<-c(-Inf,xi,Inf)
    k2<-length(xi2)
	sel<-apply(cbind(xi2[-k2],xi2[-1]),1,function(bounds) x>bounds[1]&x<=bounds[2])
	matrix(b,ncol=k+1,nrow=length(x),byrow=TRUE)[sel]
}

rPercbas<-function(n,xi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
	if (length(n)>1) n<-length(n)
	if (length(xi)!=length(F)) stop("F and xi should be of equal length")
	if (length(xi)<2) stop("The number of percentiles should be >1")
	k<-length(F)
	if (min(xi[-1]-xi[-k])<=0) stop("The elements in xi should be strictly increasing")
	if (min(F[-1]-F[-k])<=0) stop("The elements in F should be strictly increasing")
	if (F[1]!=0|F[k]!=1) stop("F[1] should be 0 and F[k] should be 1")
	qPercbas(runif(n),xi,F)
}



# The expected value, varaince and pdf of a scalar-valued X(n:r):n 
# under a percentile-based distribution
qtree.moments<-function(r,n,xi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
	newtons.binomial<-function(a,b,n) {
		pot<-seq(0,n)
		choose(n,pot)*a^(n-pot)*b^pot
	}
	npts<-length(xi)
	b<-(F[-1]-F[-npts])/(xi[-1]-xi[-npts])
	a<-F[-npts]-b*xi[-npts]
	y0=0
	x0=xi[1]
	c0<-n*choose(n-1,r-1)
	EX2<-EX<-0
	xt<-c(min(xi)-(max(xi)-min(xi))/20,min(xi))
	y1<-rep(0,2)
	for (i in 1:(npts-1)) {
		powers.A<-seq(0,r-1)
		coef.A<-newtons.binomial(a[i],b[i],r-1)
		powers.B<-seq(0,n-r)
		coef.B<-newtons.binomial(1-a[i],-b[i],n-r)
		prod<-outer(coef.A,coef.B,"*")
		powers.mat<-outer(powers.A,powers.B,"+")
		powers.vec<-unique(as.vector(powers.mat))
		polynomial<-rep(b[i],length(powers.vec))
		for (j in 1:length(powers.vec)) {
			polynomial[j]<-polynomial[j]*sum(prod[powers.mat==powers.vec[j]])
		}
		x<-seq(xi[i],xi[i+1],length=50)
		y<-c0*sapply(x,function(x) sum(polynomial*x^powers.vec))
		EX.powers<-powers.vec+2
		EX.coef<-polynomial/EX.powers
		EX<-EX+sum(EX.coef*(xi[i+1]^EX.powers-xi[i]^EX.powers))
		EX2.powers<-powers.vec+3
		EX2.coef<-polynomial/EX2.powers
		EX2<-EX2+sum(EX2.coef*(xi[i+1]^EX2.powers-xi[i]^EX2.powers))
		x<-seq(xi[i],xi[i+1],length=50)
		xt<-c(xt,x)
		y1<-c(y1,c0*b[i]*(a[i]+b[i]*x)^(r-1)*(1-a[i]-b[i]*x)^(n-r))
	}
	xt<-c(xt,max(xi),max(xi)+(max(xi)-min(xi))/20)
	y1<-c(y1,rep(0,2))
	EX<-c0*EX
	EX2<-c0*EX2
	list(mu=EX,sigma2=EX2-EX^2,x=xt,y=y1)
}

# kahden kvantiilipuun yhteisjakauma
# r1<r2
qtree.jointdens<-function(x,r1,r2,n,xi,F) {
	beta<-factorial(n)/(factorial(n-r2)*factorial(r2-r1-1)*factorial(r1-1))
	k<-length(F)
	b<-c((F[-1]-F[-k])/(xi[-1]-xi[-k]))
	a<-F[-k]-b*xi[-k]
	# fine grid over the 2-d space
	x1<-x[,1]# rep(seq(xi[1],xi[k],length=1000),1000)      # column
	x2<-x[,2]# Brep(seq(xi[1],xi[k],length=1000),each=1000) # row
	y<-rep(0,length(x1)) # density
	# go through all cells. 
	for (col in (1:(k-1))) {
		for (row in (col:(k-1))) {
			sel<-(xi[row]<x2)&(xi[row+1]>x2)&(xi[col]<x1)&(xi[col+1]>x1)
			y[sel]<-beta*b[col]*b[row]*(a[col]+b[col]*x1[sel])^(r1-1)*
					(a[row]+b[row]*x2[sel]-a[col]-b[col]*x1[sel])^(r2-r1-1)*
					(1-a[row]-b[row]*x2[sel])^{n-r2}
		} 
	}
	y[x1>=x2]<-0
	y
}

# exy: r-iimplementation on 11.4.2019
# corresponding fortran-code i n file c:/laurim/.../cov.f
# approximates the integeral of function x*y*fxy by computing for  
# each percentile interval the function mean in a regular npts*npts gridi
# and multiplying it by the area.  
qtree.exy<-function(r1,r2,n,xi,F,npts=100) {
	if (r2<=r1|r2>n|r1<1) stop("You should have n>=r2>r1>=1")
	beta<-factorial(n)/(factorial(n-r2)*factorial(r2-r1-1)*factorial(r1-1))
	k<-length(F)
	b<-c((F[-1]-F[-k])/(xi[-1]-xi[-k]))
	a<-F[-k]-b*xi[-k]
	val<-0
	# go through all cells. 
	for (col in (1:(k-1))) {
		for (row in (col:(k-1))) {
			# fine grid over the 2-d space
			xw<-diff(xi[col:(col+1)])
			yw<-diff(xi[row:(row+1)])
			stepx<-xw/npts
			stepy<-yw/npts
			x1<-rep(seq(xi[col]+stepx/2,xi[col+1],stepx),npts)      # column
			x2<-rep(seq(xi[row]+stepy/2,xi[row+1],stepy),each=npts) # row
			sel<-x1<x2
			x1<-x1[sel]
			x2<-x2[sel]
			val<-val+mean(x1*x2*beta*b[col]*b[row]*(a[col]+b[col]*x1)^(r1-1)*
							(a[row]+b[row]*x2-a[col]-b[col]*x1)^(r2-r1-1)*
							(1-a[row]-b[row]*x2)^(n-r2))*xw*yw/ifelse(col==row,2,1)
		} 
	}
	val
}



# returns tjhe pergentage values that correspond to the given expected value 
# of the quantiles, and the corresponding varaince-covaraince matrices. 
# obs is a data frame with (at least) columns r, n, plot, d.
# should be ordered by r within plot, and all observations from same plot should follow each other. 
qtree.varcov<-function(obs,xi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
	a<-matrix(unlist(apply(obs[,c("r","n")], 1,function(b) qtree.moments(b[1],b[2],xi=xi,F=F)[c(1,2)])),ncol=2,byrow=T)
	Ed<-a[,1]
	Vd<-a[,2]
	pEd<-pPercbas(Ed,xi,F=F)
	obs$Ed<-Ed
	obs$pEd<-pEd
	plots<-unique(obs$plot)
	R<-matrix(0,ncol=0,nrow=0)
	for (plot0 in plots) {
		obs0<-obs[obs$plot==plot0,]  
		Ed0<-Ed[obs$plot==plot0]
		Vd0<-Vd[obs$plot==plot0]
		R0<-diag(Vd0)
		n0<-length(Ed0)
		if (n0>1) { # compute the covarainces if needed
			for (i in 1:(n0-1)) {
				for (j in (i+1):n0) {
					R0[i,j]<-R0[j,i]<-qtree.exy(obs0$r[i],obs0$r[j],obs0$n[i],xi,F)-Ed0[i]*Ed0[j]
				} # j 
			} # i
		} # if (n0>1)
		R<-adiag(R,R0)
	}
	list(obs=obs,R=R)
}

interpolate.D<-function(D,ppi,F=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1)) {
	ppi<-round(ppi,10)
	nam<-unique(c(ppi,F))
	nam<-nam[order(nam)]
	Dapu<-D
	Fapu<-F
	for (j in 1:length(ppi)) {
		i<-max(c(1,length(Fapu[ppi[j]>=Fapu])))
#                print(c(ppi[j],Fapu[i]))
		if (ppi[j]!=Fapu[i]) {
			cov<-Dapu[i,]+(ppi[j]-Fapu[i])/(Fapu[i+1]-Fapu[i])*(Dapu[i+1,]-Dapu[i,])
			vapu<-diag(Dapu)
			var<-vapu[i]+(ppi[j]-Fapu[i])/(Fapu[i+1]-Fapu[i])*(vapu[i+1]-vapu[i])
			vcvek<-c(cov[1:i],var,cov[-c(1:i)])
			Dapu<-rbind(Dapu[1:i,],cov,Dapu[-c(1:i),])
			Dapu<-cbind(Dapu[,1:i],vcvek,Dapu[,-c(1:i)])
			Fapu<-c(Fapu[1:i],ppi[j],Fapu[-c(1:i)])
		}
	}
	D1<-matrix(nrow=length(ppi),ncol=length(ppi))
	D2<-matrix(ncol=length(ppi),nrow=length(F))
	for (i in 1:length(ppi)) {
		for (j in 1:length(ppi)) {
			D1[j,i]<-Dapu[nam==ppi[j],nam==ppi[i]]
		}
		for (j in 1:length(F)) {
			D2[j,i]<-Dapu[nam==F[j],nam==ppi[i]]
		}
	}
	list(D=Dapu,F=Fapu,D1=D1,D2=D2)
}
		 
