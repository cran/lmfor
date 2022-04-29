#### UTILITIES ####

## Transform polar coordinates to cartesian coordinates

polar_to_cart <- function(X){
  if(!is.matrix(X)) X<-matrix(X, ncol=2)
  x <- X[,1]*cos(X[,2])
  y <- X[,1]*sin(X[,2])
  
  return(cbind(x,y))
  
}

## Transform cartesian coordinates to polar coordinates

cart_to_polar <- function(X){
  if(!is.matrix(X)) X<-matrix(X, ncol=2)
  r <- sqrt(X[,1]^2 + X[,2]^2)
  phi <- atan2(X[,2],X[,1])
  
  return(cbind(r, phi))
  
}

## Calculate coordinates related to the polygons that form the area behind  discs/trees that are not visible from origin

triangle_coords <- function(X, plot.radius, polar=TRUE){
  R<-plot.radius
  if(polar){
    Y<-X[,1:2]
  }else{
    Y <- cart_to_polar(X[,1:2])
  }
  if(!is.matrix(Y)) Y<-matrix(Y, ncol=2)
  alpha <- asin(X[,3]/(2*Y[,1]))
  
  cathetus <- sqrt(Y[,1]^2 - (X[,3]/2)^2)
  
  smalltri_right <- polar_to_cart(cbind(cathetus, Y[,2]-alpha))
  smalltri_left  <- polar_to_cart(cbind(cathetus, Y[,2]+alpha))
  
  bigtri_right <- polar_to_cart(cbind(2*R, Y[,2]-alpha))
  bigtri_left  <- polar_to_cart(cbind(2*R, Y[,2]+alpha))
  
  return(cbind(smalltri_right, smalltri_left, bigtri_right, bigtri_left))
  
}



## Generate nonvisible areas produced by trees as owin objects

shades <- function(X, plot.radius, polar=TRUE){
  R<-plot.radius
  if(!is.matrix(X)) X<-matrix(X, ncol=3)
  Y <- triangle_coords(X=X, plot.radius=R, polar=polar)
  shadelist<-list()
  if(polar) X<-cbind(polar_to_cart(X[,1:2]), X[,3])
  
  for(i in 1:nrow(X)){
    
    a<-owin(poly=rbind(Y[i,1:2], Y[i,5:6], Y[i,7:8], Y[i,3:4]))
    
    a<-union.owin(a, disc(radius=X[i,3]/2, centre=c(X[i,1],X[i,2]), npoly=128))
    
    shadelist[[i]] <- a
  }
  
  return(shadelist)
}

### Helper function for preprocessing of data. Orders the rows of the data based on the closeness
### of the stem disc to the plot centre.

ordering_cps <- function(data, polar=TRUE){
  if(!is.matrix(data)){
    stop("data must be a matrix.\n")
  }
  if(nrow(data)<2){
    return(data)
  }
  
  if(polar){
    ordering<-order(data[,1]-data[,3]/2)
  }else{
    X<-matrix(cart_to_polar(data[,1:2]), ncol=2)
    ordering<-order(X[,1]-data[,3]/2)
  }
  
  
  
  data[ordering,]
  
}


#### END UTILITIES ####


## A function that categorizes trees as either detected or not detected based on a geometric visibility based detection condition

visibility_thinning_cps<-function(data, plot.radius, alpha=0, polar=TRUE, npoly=1024, delta=NULL){
  
  if(!is.matrix(data)){
    stop("data must be a matrix.\n")
    
  }else{
    if(ncol(data)!=3){
      stop("data must have 3 columns: either x, y, d or r, phi, d.\n")
      
    }
  }
  
  
  if(is.null(alpha)){
    stop("alpha should be a single number.\n")
    
  }
  if(is.null(plot.radius)){
    stop("plot.radius should be a single number.\n")
    
  }
  
  if(is.list(alpha) | is.na(alpha) |  !is.numeric(alpha)){
    stop("alpha should be a single number.\n")
    
  }
  if(is.list(plot.radius) | is.na(plot.radius) | !is.numeric(plot.radius)){
    stop("plot.radius should be a single number.\n")
    
  }
  
  
  if(length(alpha)>1){
    stop("alpha should be a single number.\n")
    
  }
  if(length(plot.radius)>1){
    stop("plot.radius should be a single number.\n")
    
  }
  
  
  R<-plot.radius
  buffer<-alpha
  
  ## Order the data based on the distance to the closest point in the stem disc
  if(polar){
    cn<-c("phi", "r", "d", "detected")
    ordering<-order(data[,1]-data[,3]/2)
  }else{
    cn<-c("x", "y", "d", "detected")
    X<-matrix(cart_to_polar(data[,1:2]), ncol=2)
    ordering<-order(X[,1]-data[,3]/2)
  }
  
  
  
  data <- matrix(data[ordering,], ncol=3)
  
  if(polar){
    X<-matrix(data[,1:2], ncol=2)
  }else{
    X<-matrix(cart_to_polar(data[,1:2]), ncol=2)
  }
  
  if(X[1,1]-data[1,3]/2 <=0){
    stop("The centre point is of the plot is covered by a stem. Visibility calculation is not possible. \n")
  }
  
  ## Produce the nonvisible areas. Categorize trees that are outside plot as not detected.
  
  n<-nrow(data)
  is.it.detected<-rep(TRUE, n)
  is.it.inside<-X[,1]<=R
  is.it.detected[!is.it.inside]<-FALSE
  
  
  shadelist<-shades(data, R, polar=polar)
  X<-matrix(polar_to_cart(data[,1:2]), ncol=2)
  
  ## Produce the polygonal approximation of the sampling window
  if(is.null(delta)){
    w.of.i <- disc(radius=R, centre=c(0,0), npoly=npoly)
  }else{
    w.of.i <- disc(radius=R, centre=c(0,0), delta=delta)
  }
  
  if(n==1){
    
    if(is.it.inside){
      treelist<-cbind(data, 1)
      colnames(treelist)<-cn
      return(treelist)
    }else{
      warning("No trees are located in the plot.\n")
      treelist<-cbind(data, 0)
      colnames(treelist)<-cn
      return(treelist)
    }
    
    
    return(cbind(data, 1))
    
  }else{
    
    a <- shadelist[[1]]
    
    ## Centre point visibility detection condition
    if(buffer==0){
      
      for(i in 2:n){
        if(is.it.inside[i]){
          b<-a
          b <- intersect.owin(b, w.of.i, fatal=F)
          if(!is.null(b)){
            if(inside.owin(x=X[i,1],y=X[i,2],w=b)){ 
              is.it.detected[i] <- FALSE
            }  
          }
        }
        a <- union.owin(a,shadelist[[i]])
      }
      
    }else{
      ## This part handles the any visibility detection condition case, and all the other erosion based cases
      if(buffer<0){
        for(i in 2:n){
          if(is.it.inside[i]){
            b <- erosion(a, abs(buffer)*data[i,3]/2)
            b <- intersect.owin(b, w.of.i, fatal=F)
            if(!is.null(b)){
              if(inside.owin(x=X[i,1],y=X[i,2],w=b)){ 
                is.it.detected[i] <- FALSE
                
              }
            }
          }
          a <- union.owin(a,shadelist[[i]])
        }
      }else{
        ## This part handles the full visibility detection condition case, and all the other dilation based cases
        for(i in 2:n){
          if(is.it.inside[i]){
            b <- dilation(a, buffer*data[i,3]/2)
            b <- intersect.owin(b, w.of.i, fatal=F)
            if(!is.null(b)){
              if(inside.owin(x=X[i,1],y=X[i,2],w=b)){ 
                is.it.detected[i] <- FALSE
                
              }
            }
          }
          a <- union.owin(a,shadelist[[i]])
        }
      }
    }
    
    
  }
  
  if(sum(is.it.inside)==0){
    warning("No trees are located in the plot.\n")
  }
  
  treelist<-cbind(data, ifelse(is.it.detected, 1, 0))
  colnames(treelist)<-cn
  
  
  treelist

}


## A function that produces detection probabilities or detectabilities in a circular sample plot

detectability_cps<-function(data, plot.radius, alpha=0, polar=TRUE, npoly=1024, delta=NULL){
  
  if(!is.matrix(data)){
    stop("data must be a matrix.\n")
    
  }else{
    if(ncol(data)!=4){
      stop("data must have 4 columns: either x, y, d, and detected, or r, phi, d and detected.\n")
      
    }
  }
  
  
  if(is.null(alpha)){
    stop("alpha should be a single number.\n")
    
  }
  if(is.null(plot.radius)){
    stop("plot.radius should be a single number.\n")
    
  }
  
  if(is.list(alpha) | is.na(alpha) |  !is.numeric(alpha)){
    stop("alpha should be a single number.\n")
    
  }
  if(is.list(plot.radius) | is.na(plot.radius) | !is.numeric(plot.radius)){
    stop("plot.radius should be a single number.\n")
    
  }
  
  
  if(length(alpha)>1){
    stop("alpha should be a single number.\n")
    
  }
  if(length(plot.radius)>1){
    stop("plot.radius should be a single number.\n")
    
  }
  
  R<-plot.radius
  buffer<-alpha
  
  ## Order the data based on the distance to the closest point in the stem disc
  if(polar){
    cn<-c("phi", "r", "d", "detected", "detectability")
    ordering<-order(data[,1]-data[,3]/2)
  }else{
    cn<-c("x", "y", "d", "detected", "detectability")
    X<-matrix(cart_to_polar(data[,1:2]), ncol=2)
    ordering<-order(X[,1]-data[,3]/2)
  }
  
  data <- matrix(data[ordering,], ncol=4)
  
  
  n<-nrow(data)
  
  detectability <- rep(1, n)
  
  
  if(polar){
    X<-matrix(data[,1:2], ncol=2)
  }else{
    X<-matrix(cart_to_polar(data[,1:2]), ncol=2)
  }
  if(X[1,1]-data[1,3]/2 <=0){
    stop("The centre point of the plot is covered by a stem. Visibility calculation is not possible. \n")
  }
  
  ## Produce the nonvisible areas. Categorize trees that are outside plot as not detected.
  shadelist<-shades(cbind(X, data[,3]), R, polar=T)
  is.it.inside<-X[,1]<=R
  detectability[!is.it.inside] <- 0
  detectability[data[,4]==0] <- 0
  
  ## Produce the polygonal approximation of the sampling window
  if(is.null(delta)){
    w.of.i <- disc(radius=R, centre=c(0,0), npoly=npoly)
  }else{
    w.of.i <- disc(radius=R, centre=c(0,0), delta=delta)
  }
  
  
  if(n==1){
    
    if(is.it.inside){
      
      treelist<-cbind(data, 1)
      colnames(treelist)<-cn
      return(treelist)
    }else{
      warning("Data includes trees that have been marked as detected but are located outside the plot boundary. The trees in question have been excluded from the estimates, detection probabilities not calculated.\n")
      treelist<-cbind(data, 0)
      colnames(treelist)<-cn
      return(treelist)
    }
    
    
  }else{
    
    a <- shadelist[[1]]
    
    ## Centre point visibility detection condition
    if(buffer==0){
      
      for(i in 2:n){
        if(data[i,4]==1 & is.it.inside[i]){
          b<-a
          b <- intersect.owin(b, w.of.i, fatal=F)
          
          c.of.i <- edges(disc(radius=data[i,1], centre=c(0,0), npoly=1024))
          detectability[i] <- 1 - sum(lengths_psp(c.of.i[b]))/(2*pi*X[i,1])
        }
        
        a <- union.owin(a,shadelist[[i]])
      }
      
    }else{
      ## This part handles the any visibility detection condition case, and all the other erosion based cases
      if(buffer<0){
        for(i in 2:n){
          if(data[i,4]==1 & is.it.inside[i]){
            b <- erosion(a, abs(buffer)*data[i,3]/2)
            b <- intersect.owin(b, w.of.i, fatal=F)
            
            
            if(is.null(delta)){
              c.of.i <- edges(disc(radius=X[i,1], centre=c(0,0), npoly=npoly))
            }else{
              c.of.i <- edges(disc(radius=X[i,1], centre=c(0,0), delta=delta))
            }
            
            detectability[i] <- 1 - sum(lengths_psp(c.of.i[b]))/(2*pi*X[i,1])
            
            
          }
          a <- union.owin(a,shadelist[[i]])
        }
      }else{
        ## This part handles the full visibility detection condition case, and all the other dilation based cases
        for(i in 2:n){
          if(data[i,4]==1 & is.it.inside[i]){
            b <- dilation(a, buffer*data[i,3]/2)
            b <- intersect.owin(b, w.of.i, fatal=F)
            
            if(is.null(delta)){
              c.of.i <- edges(disc(radius=X[i,1], centre=c(0,0), npoly=npoly))
            }else{
              c.of.i <- edges(disc(radius=X[i,1], centre=c(0,0), delta=delta))
            }
            
            detectability[i] <- 1 - sum(lengths_psp(c.of.i[b]))/(2*pi*X[i,1])
            
            
          }
          a <- union.owin(a,shadelist[[i]])
        }
      }
    }
    
    
  }
  
  
  
  treelist<-cbind(data, detectability)
  colnames(treelist)<-cn
  
  if(sum(data[,4]==1 & !is.it.inside)>0){
    warning("Data includes trees that have been marked as detected but are located outside the plot boundary. The trees in question have been excluded from the estimates, detection probabilities not calculated.\n")
  }
  
  
  treelist
  
}



## A function that produces Horvitz-Thompson-like estimates in a circular sample plot

HTest_cps<-function(data, total=TRUE, confidence.level=0.95){
  
  if(!is.matrix(data)){
    stop("data must be a matrix.\n")
  
  }else{
    if(ncol(data)<2){
      stop("data must have at least 2 columns: first column contains detectabilities, the others marks over which estimation is done.\n")
      
    }
  }
  
  ## Use only trees that have a detectability
  
  data<-matrix(data[data[,1]>0,], ncol=ncol(data))
  n<-nrow(data)
  
  if(n==0){
    warning("Data does not seem to contain detected trees.\n")
    return(0)
  }
  
  if(n==1){
    warning("Data contains only 1 detected tree.\n")
    return(cbind(data[,2:ncol(data)], rep(0, ncol(data)-1), data[,2:ncol(data)], data[,2:ncol(data)]))
  }
 
  ###calculate the right quantile to use based on the confidence level of the interval
  cl<- 1 - (1-confidence.level)/2
  
  #Produce a result matrix, where column 1 = point estimate, col 2 = variance estimate, col3 & col4 are the end points of confidence interval
  results<-matrix(0, ncol=4, nrow=ncol(data)-1)
  
  if(total){
    
    for(i in 2:ncol(data)){
      results[i-1, 1]  <- sum(data[,i]/data[,1])
      results[i-1, 2]  <- sum((1/(data[,1]^2) - 1/data[,1])*data[,i]^2)
      
      if(n<50){
        results[i-1, 3] <-  results[i-1, 1] - qt(cl, n-1)*sqrt(results[i-1, 2])
        results[i-1, 4] <-  results[i-1, 1] + qt(cl, n-1)*sqrt(results[i-1, 2])
      }
      else{
        results[i-1, 3] <-  results[i-1, 1] - qnorm(cl)*sqrt(results[i-1, 2])
        results[i-1, 4] <-  results[i-1, 1] + qnorm(cl)*sqrt(results[i-1, 2])
      }
      
    }
    
    
  }else{
    nest<-sum(1/data[,1])
    
    for(i in 2:ncol(data)){
      results[i-1, 1]  <- sum(data[,i]/data[,1])/nest
      results[i-1, 2]  <- (1/nest^2)*sum((1/(data[,1]^2) - 1/data[,1])*(data[,i]-results[i-1, 1])^2)
      
      if(n<50){
        results[i-1, 3] <-  results[i-1, 1] - qt(cl, n-1)*sqrt(results[i-1, 2])
        results[i-1, 4] <-  results[i-1, 1] + qt(cl, n-1)*sqrt(results[i-1, 2])
      }
      else{
        results[i-1, 3] <-  results[i-1, 1] - qnorm(cl)*sqrt(results[i-1, 2])
        results[i-1, 4] <-  results[i-1, 1] + qnorm(cl)*sqrt(results[i-1, 2])
      }
    } 
  }
  colnames(results)<-c("estimate", "estimated variance", "conf.low", "conf.high")
  if(!is.null(colnames(data[,-1]))) rownames(results) <- colnames(data[,-1])
  
  results
  
}

