#Scenario 2: spatial pattern driven by cell number(false positive)
#ave: mean of negative binomial distribution across all spots.
#n:index spot's mean of poisson distribution modeling cell number 
kernalcn<-function(y1,y2,x,k,r,n,l,ave,base,sz){
  set.seed(103)
  aresult<-matrix(rep(0,(l^2)),nrow=l,ncol=l)
  result<-matrix(rep(0,(2*x-1)^2),nrow=2*x-1,ncol=2*x-1)
  tresult<-matrix(rnbinom(l^2,mu=ave,size=sz)*(rpois(l^2,base)+1),nrow=l,ncol=l)
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
#randam generate cell number following poission distribution with mean  
  generate2<-function(s){
    numb<-rpois(1,lambda=s)+1
    g<-sum(rnbinom(numb,mu=ave,size=sz))
    return(list(numb,g))
  }
  cresult<-sapply(result,generate2)
  #dresult is sum counts of each spot in pattern part
  dresult<-matrix(as.numeric(cresult[2,]),nrow=2*x-1,ncol=2*x-1,byrow=T)
  #nresult is cell number of each spot
  nresult<-matrix(as.numeric(cresult[1,]),nrow=2*x-1,ncol=2*x-1,byrow=T)

  # embedd pattern into background  
  aresult[y1,y2]<-nresult[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    aresult[y1-i,y2]<-nresult[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult[y1-i,y2-j]<-nresult[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult[y1-i,y2+t]<-nresult[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    aresult[y1+s,y2]<-nresult[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult[y1+s,y2-j]<-nresult[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult[y1+s,y2+t]<-nresult[x+s,x+t]
    }
  }
  
  #spot level pattern
  tresult[y1,y2]<-dresult[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    tresult[y1-i,y2]<-dresult[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      tresult[y1-i,y2-j]<-dresult[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      tresult[y1-i,y2+t]<-dresult[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    tresult[y1+s,y2]<-dresult[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      tresult[y1+s,y2-j]<-dresult[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      tresult[y1+s,y2+t]<-dresult[x+s,x+t]
    }
  }
  #Pattern of Average gene expression  
  avreuslt<-tresult/aresult
  
  graph1<-filled.contour(x = 1:nrow(aresult),y = 1:ncol(aresult),
                         z = aresult, color.palette = myPalette,
                         plot.title = title(main = "Distribution of Cell Number",
                                            xlab = "x-coordinate",ylab = "y-coordinate"),
                         plot.axes = {axis(1, seq(1, ncol(aresult), by = 5))
                           axis(2, seq(1, nrow(aresult), by = 5))},
                         key.title = title(main="Cell\n(number)"),
                         key.axes = axis(4, seq(min(aresult), max(aresult), by = 10))
                         ,)
  graph2<-filled.contour(x = 1:nrow(tresult),y = 1:ncol(tresult),
                         z = tresult, color.palette = myPalette,
                         plot.title = title(main = "Pattern of Driven by Cell Number",
                                            xlab = "x-coordinate",ylab = "y-coordinate"),
                         plot.axes = {axis(1, seq(1, ncol(tresult), by = 5))
                           axis(2, seq(1, nrow(tresult), by = 5))},
                         key.title = title(main="Gene\n(counts)"),
                         key.axes = axis(4, seq(min(tresult), max(tresult), by = 1000))
                         ,)
  graph3<-filled.contour(x = 1:nrow(avreuslt),y = 1:ncol(avreuslt),
                         z = avreuslt, color.palette = myPalette,
                         plot.title = title(main = "Pattern of Average Gene Expression",
                                            xlab = "x-coordinate",ylab = "y-coordinate"),
                         plot.axes = {axis(1, seq(1, ncol(avreuslt), by = 5))
                           axis(2, seq(1, nrow(avreuslt), by = 5))},
                         key.title = title(main="Gene\n(counts)"),
                         key.axes = axis(4, seq(min(avreuslt), max(avreuslt), by = 3))
                         ,)
  finalresult<-list(graph1,graph2,graph3)
  return(finalresult)
}
kernal5<-kernalcn(30,25,40,2,1.1,100,50,90,5,10)
