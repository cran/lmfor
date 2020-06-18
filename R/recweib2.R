# TODO: Add comment
# 
# Author: Lauri Mehtatalo 20.2.2020
###############################################################################

# R-script for parameter recovery for the basal-area weighted 2-parameter Weibull function assumed for diameter distribution
# Assume we know one of the following mean characteristics:
# A: D = arithmetic mean diameter i.e. the first moment
# B: DG = basal area weighted mean diameter i.e. mean of the weighted second order Weibull distribution 
# C: DM = arithmetic median diameter i.e. 0.5 quantile
# D: DGM = basal area median diameter i.e. 0.5 quantile of the weighted second order Weibull distribution
# In addition N and G  are known

# Basal area (G) is given as m^2/ha
# The number of stems (N) is given as 1/ha
# The mean diameter is given in cm.


# Functions returning the Weibull scale parameter for given shape parameter and mean diameter, 
# for the four altermnative mean dianmeters.

# A: Unweighted mean
scaleDMean2<-function(D,shape) {
	D*gamma(1-2/shape)/gamma(1-1/shape) 
}

# B Basal-area weighted mean    
scaleDGMean2<-function(D,shape) {
	D/gamma(1/shape+1) 
}

# C: Unweighted median
scaleDMed2<-function(D,shape) {
	D/(qgamma(0.5,1-2/shape))^{1/shape}      
}

# D: Basal-area weighted median
scaleDGMed2<-function(D,shape) {
	D/(log(2)^(1/shape))   
}  

# The function to be set to zero: 
# 40000*G/(pi*N)-scale^2/gamma(1-2/shape).
# Dtype should be any of letters "A", "B", "C" and "D",
# And specifies which type of mean or median diameter was given 
# in argument "D".

func.recweib2<-function(lshape,G,N,D,Dtype, trace=FALSE) {  
	shape<-exp(lshape)+2.01
	if (Dtype=="A") {
		scale<-scaleDMean2(D,shape)
	} else if (Dtype=="B") {
		scale<-scaleDGMean2(D,shape)
	} else if (Dtype=="C") {
		scale<-scaleDMed2(D,shape)
	} else if (Dtype=="D") {
		scale<-scaleDGMed2(D,shape)
	} else stop("Dtype should be any of 'A', 'B', 'C', 'D'")
	
	sol<-40000*G/(pi*N)-scale^2/gamma(1-2/shape) 
	attributes(sol)$scale<-scale
	if (trace) cat(shape,attributes(sol)$scale,sol,"\n")
	sol
}

# The function for recovery of shape and scale.
# G is Basal area, m^2/ha
# N Number of stems, 1/ha
# Weighted or unweighted D is mean or median diameter, cm
# Dtype the type of the given mean diameter, given as "A", "B", "C" and "D".
# A: The arithmetic mean 
# B: The basal-area weighted mean
# C: The median 
# D: The median weighted by basal area
# trace: if TRUE, then intermediate output during estimation is printed on the screen.
# init: the initial guess for the shape parameter. 
# If not given, value shape=4 is used. 
# 
# The function returns a list with the following elements
# shape: the shape parameter (c) in the solution
# scale: the scale parameter (b) in the solution 
# G, N, D, Dtype: the values of the input arguments
# val: the value of the equation shape(G,N,D)-shape at the solution
# 
recweib2<-function(G,N,D,Dtype,init=NA,trace=FALSE) {
	lshape<-NRnum(init=log(init-2.01), 
			fnlist=list(function(lshape) func.recweib2(lshape,G=G,N=N,D=D,Dtype=Dtype,trace=trace)),
			crit=10)
	shape=exp(lshape$par)+2.01
	if (Dtype=="A") {
		scale<-scaleDMean2(D,shape)
	} else if (Dtype=="B") {
		scale<-scaleDGMean2(D,shape)
	} else if (Dtype=="C") {
		scale<-scaleDMed2(D,shape)
	} else if (Dtype=="D") {
		scale<-scaleDGMed2(D,shape)
	} 
	
	list(shape=shape, scale=scale, G=G, N=N, D=D, Dtype=Dtype, val=lshape$value)
}

# TODO: Add comment
# 
# Author: Lauri Mehtatalo and Jouni Siipilehto 27.1.2014
###############################################################################

# Siipilehto & Meht2talo 2013. "Parameter recovery vs. parameter prediction for the Weibull distribution
# validated for Scots pine stands in Finland". Silva Fennica 47(4), article id 1057; 

# R-script for parameter recovery for the unweighted 2-parameter Weibull function assumed for diameter distribution
# Assume we know one of the following mean characteristics:
# A: D = arithmetic mean diameter i.e. the first moment
# B: DG = basal area weighted mean diameter i.e. mean of the weighted second order Weibull distribution 
# C: DM = arithmetic median diameter i.e. 0.5 quantile
# D: DGM = basal area median diameter i.e. 0.5 quantile of the weighted second order Weibull distribution
# The second moment (Dq^2) is calculated from the stem number (N_ha) and basal area per hectare (G_ha)

# Basal area (G) is given as m^2/ha
# The number of stems (N) is given as 1/ha
# The mean diameter is given in cm.


# Functions returning the Weibull scale parameter for given shape and mean diameter, 
# for the four altermnative mean dianmeters.
# Functions DXXX return the first derivative of these functions, 
# but they are not currently in use and may have mistakes. 

# A: Unweighted mean
scaleDMean1<-function(D,shape) {
	D/gamma(1/shape+1) 
}

# B Basal-area weighted mean    
scaleDGMean1<-function(D,shape) {
	D*gamma(2/shape+1)/gamma(3/shape+1) 
}


# C: Unweighted median
scaleDMed1<-function(D,shape) {
	D/(log(2)^(1/shape))   
}  


# D: Basal-area weighted median
scaleDGMed1<-function(D,shape) {
	D/(qgamma(0.5,2/shape+1))^{1/shape}      
}

# The function to be set to zero: 
# Returns the difference of the given shape parameter and 
# the value corresponding to given stand chaacteristics.
# Dtype should be any of letters "A", "B", "C" and "D",
# And specifies which type of mean or median diameter was given 
# in argument "D".

func.recweib1<-function(lshape,G,N,D,Dtype, trace=FALSE) {  
	shape<-exp(lshape)+0.01
	if (Dtype=="A") {
		scale<-scaleDMean1(D,shape)
	} else if (Dtype=="B") {
		scale<-scaleDGMean1(D,shape)
	} else if (Dtype=="C") {
		scale<-scaleDMed1(D,shape)
	} else if (Dtype=="D") {
		scale<-scaleDGMed1(D,shape)
	} else stop("Dtype should be any of 'A', 'B', 'C', 'D'")
	
	sol<-40000*G/(pi*N)-scale^2*gamma(2/shape+1) 
	attributes(sol)$scale<-scale
	if (trace) cat(shape,attributes(sol)$scale,sol,"\n")
	sol
}



# The function for recovery of shape and scale.
# G is Basal area, m^2/ha
# N Number of stems, 1/ha
# Weighted or unweighted D is mean or median diameter, cm
# Dtype the type of the given mean diameter, given as "A", "B", "C" and "D".
# A: The arithmetic mean 
# B: The basal-area weighted mean
# C: The median 
# D: The median weighted by basal area
# trace: if TRUE, then intermediate outbut during estimation is printed on the screen.
# init: the initial guess for the shape parameter. 
# If not given, the initial guesses are based on the PPM equations given 
# in the Appendix of Siipilehto and Mehtatalo (2013).
# 
# The recovery is based on the solution of the equation 
# shape(G,N,D)-shape = 0
# Where shape(G,N,D) expresses the shape parameter of The weibull distribution
# for the given values of the stand characteristics and shape is the estimated
# value of the shape parameter.
# 
# The function returns a list with the following elements
# shape: the shape parameter (c) in the solution
# scale: the scale parameter (b) in the solution 
# G, N, D, Dtype: the values of the input arguments
# val: the value of the equation shape(G,N,D)-shape at the solution
# 
recweib1<-function(G,N,D,Dtype,init=NA,trace=FALSE) {
	lshape<-NRnum(init=log(init-0.01), 
			fnlist=list(function(lshape) func.recweib1(lshape,G=G,N=N,D=D,Dtype=Dtype,trace=trace)),
			crit=10)
	shape=exp(lshape$par)+0.01
	if (Dtype=="A") {
		scale<-scaleDMean1(D,shape)
	} else if (Dtype=="B") {
		scale<-scaleDGMean1(D,shape)
	} else if (Dtype=="C") {
		scale<-scaleDMed1(D,shape)
	} else if (Dtype=="D") {
		scale<-scaleDGMed1(D,shape)
	} 
	
	list(shape=exp(lshape$par)+0.01, scale=scale, G=G, N=N, D=D, Dtype=Dtype, val=lshape$value)
}

recweib<-function(G,N,D,Dtype,init=NA,trace=FALSE,weight=0) {
         if (!sum(Dtype==c("A","B","C","D"))) stop("Dtype should be any of 'A', 'B', 'C', 'D'")
         if (!sum(weight==c(0,2))) stop("weight should be either 0 or 2")
         if (weight==0) {
            DQM<-sqrt(40000*G/(pi*N)) 
            if (is.na(init)) {
               if (Dtype=="A") {
                  init<-1.0/log(DQM^4/D^4 + 0.1)
                  } else if (Dtype=="B") {
                  init<-1/log(D^2/DQM^2 + 0.05)
                  } else if (Dtype=="C") {
                  init<-1/log(DQM^4/D^4 + 0.1)
                  } else if (Dtype=="D") {
                  init<-1/log(D^2/DQM^2 + 0.05)
                  } 
               } # (is.na(init))
               recweib1(G,N,D,Dtype,init,trace)
            } else if (weight==2) {
              if (is.na(init)) init<-4 
              recweib2(G,N,D,Dtype,init,trace)
              }
         }


## A function that computes different stand characteristics using Wibull distribution with parameters shape and scale
## Just to check the results
#		 chars<-function(shape,scale,G) {
#			 max<-qweibull(1-1e-10,shape,scale)
#			 # G/N
#			 GperN<-pi/40000/integrate(function(x) x^(-2)*dweibull(x,shape,scale),0,max)$value
#			 DGmed<-qweibull(0.5,shape,scale)
#			 DGmean<-integrate(function(x) x*dweibull(x,shape,scale),0,max)$value
#			 Dmean<-scale*gamma(1-1/shape)/gamma(1-2/shape)
#			 Dmed<-scale*(qgamma(0.5,1-2/shape))^{1/shape}      
#			 list(gperN=GperN,N=G/GperN,DGmed=DGmed,DGmean=DGmean,Dmean=Dmean,Dmed=Dmed)
#		 }
#		 
# Demonstration with one example stand.
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

#		 
## Example 2. Even aged stand in Finland (see Siipilehto & Mehtatalo, Fig 2):
#		 G_ha<-9.6
#		 N_ha<-949
#		 D<-11.0
#		 DG<-12.3
#		 DM<-11.1
#		 DGM<-12.4  
#		 recweib2(G_ha,N_ha,D,"A") # 4.465673 12.05919
#		 recweib2(G_ha,N_ha,DG,"B") # 4.463991 12.05912
#		 recweib2(G_ha,N_ha,DM,"C")  # 4.410773 12.05949
#		 recweib2(G_ha,N_ha,DGM,"D") # 4.448272 12.05924
#		 
#		 
## Example 3. Assumed peaked even aged stand (see Siipilehto & Mehtatalo, Fig 1):
#		 G_ha<-10.0
#		 N_ha<-1300
#		 D<-9.89
#		 DG<-10.0
#		 DM<-9.89
#		 DGM<-10.0  
#		 a<-recweib2(G_ha,N_ha,D,"A")  #  34.542 10.04978
#		 b<-recweib2(G_ha,N_ha,DG,"B") # 14.23261 10.22781
#		 c<-recweib2(G_ha,N_ha,DM,"C") # 6.708882 10.44448
#		 d<-recweib2(G_ha,N_ha,DGM,"D") # 24.45228 10.10607
#		 
#		 unlist(chars(a$shape,a$scale,G_ha))
#		 unlist(chars(b$shape,b$scale,G_ha))
#		 unlist(chars(c$shape,c$scale,G_ha))
#		 unlist(chars(d$shape,d$scale,G_ha))
#		 