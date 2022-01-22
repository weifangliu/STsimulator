#start
library(RColorBrewer)
library(gplots)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#Hotspot：y1,y2:coordinates of index spot
#2x-1:length of rectangular with pattern, x should be odd number. x-1 need to be divide by k.
#k:thinkness of every layer
#r:change rate
#n: peak size parameter in alternative parameterization of NB. Variance=mu+mu^2/size
#n1: mean parameter in NB
#l:length of background rect.
kernal1<-function(y1,y2,x,k,r,n1,n,l){
  set.seed(104)
  aresult<-matrix(rep(0,(l^2)),nrow=l,ncol=l)
  result<-matrix(rep(0,(2*x-1)^2),nrow=2*x-1,ncol=2*x-1)
  result[x,x]<-n
  # x-1 has to be divided by k
  #four corner on every k "layer"
  for (t in (0:((x-1)/k-1))) {
    result[x-(t+1)*k,x+(t+1)*k]<-(1/r)*result[x-t*k,x+t*k]
    result[x+(t+1)*k,x+(t+1)*k]<-(1/r)*result[x+t*k,x+t*k]
    result[x-(t+1)*k,x-(t+1)*k]<-(1/r)*result[x-t*k,x-t*k]
    result[x+(t+1)*k,x-(t+1)*k]<-(1/r)*result[x+t*k,x-t*k]
    #four corner on every layer between k-1 and k    
    for (i in (0:(k-1))) {
      result[x-(t+1)*k+i,x+(t+1)*k-i]<-result[x-(t+1)*k,x+(t+1)*k]
      result[x+(t+1)*k-i,x+(t+1)*k-i]<-result[x+(t+1)*k,x+(t+1)*k]
      result[x-(t+1)*k+i,x-(t+1)*k+i]<-result[x-(t+1)*k,x-(t+1)*k]
      result[x+(t+1)*k-i,x-(t+1)*k+i]<-result[x+(t+1)*k,x-(t+1)*k]
      #every line aligned to the corner (clockwisely)  
      for (j in (1:(2*((t+1)*k-i)-1))) {
        result[x-(t+1)*k+i+j,x+(t+1)*k-i] <-result[x-(t+1)*k+i,x+(t+1)*k-i]
        result[x+(t+1)*k-i,x+(t+1)*k-i-j] <-result[x+(t+1)*k-i,x+(t+1)*k-i]
        result[x-(t+1)*k+i,x-(t+1)*k+i+j] <-result[x-(t+1)*k+i,x-(t+1)*k+i]
        result[x+(t+1)*k-i-j,x-(t+1)*k+i] <-result[x+(t+1)*k-i,x-(t+1)*k+i]
      }
    }
  }
  generate<-function(n){
    g<-mean(rnbinom(5,mu=n1,size=n))
  }
  cresult<-sapply(result,generate)
  dresult<-matrix(cresult,nrow=2*x-1,ncol=2*x-1,byrow=T)
  
  aresult[y1,y2]<-dresult[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    aresult[y1-i,y2]<-dresult[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult[y1-i,y2-j]<-dresult[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult[y1-i,y2+t]<-dresult[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    aresult[y1+s,y2]<-dresult[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult[y1+s,y2-j]<-dresult[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult[y1+s,y2+t]<-dresult[x+s,x+t]
    }
  }

  finalresult<-aresult
  return(finalresult)
}
#example of kernal1
k11<-kernal1(20,30,60,2,1.1,7,250,65)



#linear gradient spot-level pattern
#y2:index column location
#2y-1: length of pattern
#r:change rate
#k:thinkness of each layer
#n1:expected mean of gene counts of spot in index column
#n:expected size of NB
#l: length of background
kernal2<-function(y2,y,r,k,n1,n,l){
  set.seed(103)
  result<-matrix(rep(0,l*((2*y)-1)),nrow=l,ncol=2*y-1)
  aresult<-matrix(rep(0,l^2),nrow=l,ncol=l)
  result[,y]<-n
  for (t in (0:((y-1)/k-1))){
    result[,y-(t+1)*k]<-(1/r)*result[,y-t*k]
    result[,y+(t+1)*k]<-(1/r)*result[,y+t*k]
    for (j in (0:(k-1))){
      result[,y-(t*k)+j]<-result[,y-t*k]
      result[,y+(t*k)-j]<-result[,y+t*k]
    }
  }
  #  10,20,1.1,1,100,2,30
  generate<-function(s){
    g<-mean(rnbinom(20,mu=n1,size=s))
  } 
  cresult<-sapply(result,generate)
  dresult<-matrix(cresult,nrow=l,ncol=(2*y)-1,byrow=F)
  aresult[,y2]<-dresult[,y]
  for (i in (0:min(c(y-1,y2-1)))){
    aresult[,y2-i]<-dresult[,y-i]
  }
  for (j in (0:min(c(l-y2,y-1)))){
    aresult[,y2+j]<-dresult[,y+j]
  }
  graph<-filled.contour(x = 1:ncol(aresult),y = 1:nrow(aresult),
                        z = t(aresult), color.palette = myPalette,
                        plot.title = title(main = "Linear Gradient Spot-level Pattern ",
                                           xlab = "x-coordinate",ylab = "y-coordinate"),
                        plot.axes = {axis(1, seq(1, ncol(aresult), by = 5))
                          axis(2, seq(1, nrow(aresult), by = 5))},
                        key.title = title(main="Gene\n(counts)"),
                        key.axes = axis(4, seq(min(aresult), max(aresult), by = 10)),
  )
  finalresult<-aresult
  return(finalresult)
}
k12<-kernal2(20,50,1.05,1,100,2,65)

#streak
kernal3<-function(y1,y2,n1,n,l,base){
  result<-matrix(rep(0,l^2),nrow=l,ncol=l)
  result[,y1:y2]<-n
  result[,0:y1]<-base
  result[,y2:l]<-base
  generate<-function(s){
    g<-rnbinom(1,mu=n1,size=s)
  }
  cresult<-sapply(result,generate)
  aresult<-matrix(cresult,nrow=l,ncol=l,byrow=F)

  finalresult<-aresult
  return(finalresult)
}
k13<-kernal3(10,20,10,100,65,50)

