#### Help function to make the calculation of window area reasonable
area_esh <- function(W){
  
  if(is.null(W)){
    return(0)
  }
  else{
    return(area.owin(W))
  }
  
}

### help function to form a union of discs without the function 'discs', which requires ppp objects,
### which require that all center points are inside some window (a restriction we do not want for minus-sampling reasons)
gg_wind <- function(treelist){
  
  if(is.matrix(treelist)){
    
    W <- disc(radius=treelist[1,3], centre=c(treelist[1, 1], treelist[1,2]), npoly=1000)
    
    if(nrow(treelist)>1){
      for(i in 2:nrow(treelist)){
        W <- union.owin(W, disc(radius=treelist[i,3], centre=c(treelist[i, 1], treelist[i,2]), npoly=1000))
      }
    }
  }else{
    treelist<-matrix(treelist,ncol=3)
    
    W <- disc(radius=treelist[1,3], centre=c(treelist[1,1], treelist[1,2]), npoly=1000)
  } 
  
  
  return(W)
}

#### Actual HT-estimating function

HTest <- function(treelist, plotwindow, alpha){
  
  # Error handling
  if(alpha < -1 | alpha > 1){
    stop("alpha must be between -1 and 1")
    
  }
  if(!is.matrix(treelist)){
    stop("treelist must be a matrix with 3 columns: x and y coordinates and the radius")
    
  }
  
  if(!is.owin(plotwindow)){
    stop("plotwindow must be an object of class owin")
    
  }
  
  #Find the trees that are inside the window, for which the detectabilities are calculated. Sort them.
  calcpoints<-treelist[inside.owin(treelist[,1], treelist[,2], plotwindow), 3]
  calcpoints<-sort(calcpoints, decreasing = T)
  n<-length(calcpoints)
  #The vector where detectabilities are stored
  detectability <- rep(0, n)
  
  # abbreviation for the window object
  a<-plotwindow
  # area of window
  plotsize<-area.owin(a)
  
  if(alpha==0){
    # Form the first W, the set that always contains the tree crown discs larger than the tree whose detectability is 
    # being calculated
    W<-gg_wind(treelist[treelist[,3]>=calcpoints[1],])
    for (i in 2:n){
      # Check if new trees can be added to W: those that are larger than the current tree but are not in W yet
      if(sum(treelist[,3]>calcpoints[i] & treelist[,3]<=calcpoints[i-1])>0){
        tmp<-gg_wind(treelist[treelist[,3]>calcpoints[i] & treelist[,3]<=calcpoints[i-1],])
        W<-union.owin(W, tmp)
      }
      # Calculate 1 - detectability
      detectability[i] <- area_esh(intersect.owin(W, a, fatal=FALSE))/plotsize
    }
    
    
  }else{
    # Erosion case
    if(alpha>0){
      # Form the first W, the set that always contains the tree crown discs larger than the tree whose detectability is 
      # being calculated
      W<-gg_wind(treelist[treelist[,3]>=calcpoints[1],])
      for (i in 2:n){
        # Check if new trees can be added to W: those that are larger than the current tree but are not in W yet
        if(sum(treelist[,3]>calcpoints[i] & treelist[,3]<=calcpoints[i-1])>0){
          tmp<-gg_wind(treelist[treelist[,3]>calcpoints[i] & treelist[,3]<=calcpoints[i-1],])
          W<-union.owin(W, tmp)
        }
        # Calculate 1 - detectability
        detectability[i] <- area_esh(intersect.owin(erosion.owin(W,alpha*calcpoints[i], shrink.frame=FALSE), a, fatal=FALSE))/plotsize
        
      }
    }else{
      # Dilation case
      # Form the first W, the set that always contains the tree crown discs larger than the tree whose detectability is 
      # being calculated
      W<-gg_wind(treelist[treelist[,3]>=calcpoints[1],])
      for (i in 2:n){
        # Check if new trees can be added to W: those that are larger than the current tree but are not in W yet
        if(sum(treelist[,3]>calcpoints[i] & treelist[,3]<=calcpoints[i-1])>0){
          tmp<-gg_wind(treelist[treelist[,3]>calcpoints[i] & treelist[,3]<=calcpoints[i-1],])
          W<-union.owin(W, tmp)
        }
        # Calculate 1 - detectability
        detectability[i] <- area_esh(intersect.owin(dilation.owin(W,abs(alpha)*calcpoints[i]), a, fatal=FALSE))/plotsize
      }
    }  
  }
  # Calculate detectability
  detectability <- 1 - detectability
  
  # Calculate the estimated number of trees
  N <- sum(1/detectability)
  # Produce a treelist
  treelist<-cbind(calcpoints, detectability)
  colnames(treelist)<-c("r", "detectability")
  
    
    return(list('N'=N, 'treelist'=treelist))
  
  
}
