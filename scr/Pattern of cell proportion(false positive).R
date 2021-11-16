#Scenario 3:Spatial pattern driven by cell proportion(false positive)
library(scatterpie)
kernalcp<-function(y1,y2,x,k,r,n1,l,ave1,ave2,base,sz){
  set.seed(103)
#background map
  aresult1<-matrix(rep(0,(l^2)),nrow=l,ncol=l)
  aresult2<-matrix(rep(0,(l^2)),nrow=l,ncol=l)
#result is pattern on cell proportion of cell type 1 and 2
  result<-matrix(rep(0,(2*x-1)^2),nrow=2*x-1,ncol=2*x-1)
  result2<-matrix(rep(0,(2*x-1)^2),nrow=2*x-1,ncol=2*x-1)
#base gene counts of spot level of cell type 1 and 2
  tresult1<-matrix(rnbinom(l^2,mu=ave1,size=sz)*(rpois(l^2,base)+1),nrow=l,ncol=l)
  tresult2<-matrix(rnbinom(l^2,mu=ave2,size=sz)*(rpois(l^2,base)+1),nrow=l,ncol=l)
#base cell number on spot level in pattern area
  numresult<-matrix(rpois((2*x-1)^2,base)+1,nrow=2*x-1,ncol=2*x-1)
#base cell number on spot level in whole area
  Numresult<-matrix(rpois(l^2,base)+1,nrow=l,ncol=l)
  
  result[x,x]<-n1
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
# cell proportion of cell type 2
  result2<-1-result

  generate1<-function(s){
    numb<-rpois(1,lambda=s)+1
    g<-sum(rnbinom(numb,mu=ave1,size=sz))
    return(list(numb,g))
  }
  generate2<-function(s){
    numb<-rpois(1,lambda=s)+1
    g<-sum(rnbinom(numb,mu=ave2,size=sz))
    return(list(numb,g))
  }

  cresult<-sapply(result*numresult,generate1)
  cresult2<-sapply(result2*numresult,generate2)
  
  #dresult is sum counts of each spot in pattern part (cell type1)
  dresult<-matrix(as.numeric(cresult[2,]),nrow=2*x-1,ncol=2*x-1,byrow=T)
  #nresult is cell number of each spot (cell type1)
  nresult<-matrix(as.numeric(cresult[1,]),nrow=2*x-1,ncol=2*x-1,byrow=T)
  
  #dresult2 is sum counts of each spot in pattern part (cell type2)
  dresult2<-matrix(as.numeric(cresult2[2,]),nrow=2*x-1,ncol=2*x-1,byrow=T)
  #nresult2 is cell number of each spot (cell type2)
  nresult2<-matrix(as.numeric(cresult2[1,]),nrow=2*x-1,ncol=2*x-1,byrow=T)
  
#pattern on cell number (cell type 1)  
  aresult1[y1,y2]<-nresult[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    aresult1[y1-i,y2]<-nresult[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult1[y1-i,y2-j]<-nresult[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult1[y1-i,y2+t]<-nresult[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    aresult1[y1+s,y2]<-nresult[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult1[y1+s,y2-j]<-nresult[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult1[y1+s,y2+t]<-nresult[x+s,x+t]
    }
  }

#pattern on cell number (cell type 2)  
  aresult2[y1,y2]<-nresult2[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    aresult2[y1-i,y2]<-nresult2[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult2[y1-i,y2-j]<-nresult2[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult2[y1-i,y2+t]<-nresult2[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    aresult2[y1+s,y2]<-nresult2[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      aresult2[y1+s,y2-j]<-nresult2[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      aresult2[y1+s,y2+t]<-nresult2[x+s,x+t]
    }
  }  
#spot level pattern (cell type 1)
  tresult1[y1,y2]<-dresult[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    tresult1[y1-i,y2]<-dresult[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      tresult1[y1-i,y2-j]<-dresult[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      tresult1[y1-i,y2+t]<-dresult[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    tresult1[y1+s,y2]<-dresult[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      tresult1[y1+s,y2-j]<-dresult[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      tresult1[y1+s,y2+t]<-dresult[x+s,x+t]
    }
  }
  
#spot level pattern (cell type 2)
  tresult2[y1,y2]<-dresult2[x,x]
  for (i in (0:min(c((y1-1),(x-1))))){
    tresult2[y1-i,y2]<-dresult2[x-i,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      tresult2[y1-i,y2-j]<-dresult2[x-i,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      tresult2[y1-i,y2+t]<-dresult2[x-i,x+t]
    }
  }
  for (s in (0:min(c((l-y1),(x-1))))){
    tresult2[y1+s,y2]<-dresult2[x+s,x]
    for (j in (0:min(c((y2-1),(x-1))))){
      tresult2[y1+s,y2-j]<-dresult2[x+s,x-j]
    }
    for (t in (0:min(c((l-y2),(x-1))))){
      tresult2[y1+s,y2+t]<-dresult2[x+s,x+t]
    }
  }
  #Pattern of Average gene expression  
  avreuslt1<-tresult1/aresult1
  avreuslt2<-tresult2/aresult2
  
  #Spot level pattern
  tresultf<-tresult1+tresult2
  avresultf<-tresultf/(aresult1+aresult2)
    
  graph1<-filled.contour(x = 1:nrow(tresultf),y = 1:ncol(tresultf),
                         z = tresultf, color.palette = myPalette,
                         plot.title = title(main = "Pattern Drive by Cell Proportion",
                                            xlab = "x-coordinate",ylab = "y-coordinate"),
                         plot.axes = {axis(1, seq(1, ncol(tresultf), by = 5))
                           axis(2, seq(1, nrow(tresultf), by = 5))},
                         key.title = title(main="Gene\n(counts)"),
                         key.axes = axis(4, seq(min(tresultf), max(tresultf), by = 200))
                         )

  # graph2<-filled.contour(x = 1:nrow(avresultf),y = 1:ncol(avresultf),
  #                        z = avresultf, color.palette = myPalette,
  #                        plot.title = title(main = "Pattern of Average Gene Expression",
  #                                           xlab = "x-coordinate",ylab = "y-coordinate"),
  #                        plot.axes = {axis(1, seq(1, ncol(avresultf), by = 5))
  #                          axis(2, seq(1, nrow(avresultf), by = 5))},
  #                        key.title = title(main="Gene\n(counts)"),
  #                        key.axes = axis(4, seq(min(avresultf), max(avresultf), by = 3))
  #                        )
  
  k1<-as.numeric(aresult1)
  k2<-as.numeric(aresult2)
 
  cellnum<-data.frame(xaxis=rep(1:l,l), yaxis=rep(1:l,each=l))
  n<-nrow(cellnum)
  cellnum$region <- factor(1:n)
  cellnum$A <- k1
  cellnum$B <- k2
  cellnum$radius<- 0.31
  colnames(cellnum)[3+(1:length(num))]<-letters[1:length(num)]
  
  p <- ggplot() + geom_scatterpie(aes(x=xaxis, y=yaxis, group=region, r=radius), data=cellnum,
                                  cols=LETTERS[1:2], color=NA) + coord_equal()
  p + geom_scatterpie_legend(cellnum$radius, x=0, y=0)
  finalresult<-list(graph1,graph2,aresult1,aresult2)
  return(finalresult)
}
library(ggplot2)
kernal6<-kernalcp(10,20,31,2,1.1,1,30,100,50,50,20)

k1<-as.numeric(kernal6[[3]])
k2<-as.numeric(kernal6[[4]])

cellnum<-data.frame(xaxis=rep(1:30,30), yaxis=rep(1:30,each=30))
n<-nrow(cellnum)
cellnum$region <- factor(1:n)
cellnum$A <- k1
cellnum$B <- k2
cellnum$radius<- 0.31
colnames(cellnum)[2+(1:length(num))]<-letters[1:length(num)]

p <- ggplot() + geom_scatterpie(aes(x=xaxis, y=yaxis, group=region, r=radius), data=cellnum,
                                cols=LETTERS[1:2], color=NA) + coord_equal()
p + geom_scatterpie_legend(cellnum$radius, x=0, y=0)

library(ggplot2)
cols = c("#90BFF9", "red")
df <- data.frame(xval=c(1:100), yval=c(1:100))
test.plot = ggplot(df, aes(x=xval, y=yval, colour=yval)) + scale_color_gradientn(colors = cols, guide = "colorbar")
col.used = ggplot_build(test.plot)$data[[1]][,1]
