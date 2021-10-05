#start
library(RColorBrewer)
library(gplots)
library(MASS)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

#Hotspot：y1,y2:coordinates of index spot
#2x-1:length of rectangular with pattern, x should be odd number. x-1 need to be divide by k.
#k:thinkness of every layer
#r:change rate
#n: peak correlation
#l:length of background rect.
#mn1,mn2:mean of biviriate normal
#v1,v2:variance of biviriate normal (assuming unequal variance)
kernalv1<-function(y1,y2,x,k,r,n,mn1,mn2,v1,v2,l){
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
#correlation pattern
  generate<-function(r){
    cov<-matrix(nrow=2,ncol=2)
    diag(cov)<-c(v1,v2)
    cov[1,2]<-r*sqrt(cov[1,1])*sqrt(cov[2,2])
    cov[2,1]<-cov[1,2]
    g<-mvrnorm(n=1,mu=as.vector(c(mn1,mn2)),Sigma=cov)
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
  finalresult<-list(aresult)
  return(finalresult)
}
#example of kernalv1
k11<-kernalv1(20,30,60,2,1.1,1,10,20,4,5,65)
