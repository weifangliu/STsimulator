#Simulation of spot-level pattern (hotspot)


#In this simulation model, gene expression counts of each spot (sum of gene counts of all cells within each spot) is modeled by negative binomial distribution. It allows to specify desired spatial hotspot-like patterns at spot level. One gene is considered in this simulation model


#' @param y1 first coordinate of index spot
#' @param y2 second coordinate of index spot
#' @param k thinkness of every layer
#' @param x length of rectangular with pattern, x should be odd number. \code{x-1} need to be divide by \code{k}
#' @param r change rate, the outer layer's mean is \code{1/r} times of that of the adjacent inner layer
#' @param n mean parameter of the negative binomial distribution of index spot
#' @param sz size parameter in alternative parametrization of negative binomial distribution
#' @param l length of background

#'
#' @return A graph with hotspot-like pattern of SVG at spot level

#' @useDynlib sost
#' @importFrom gplots RColorBewer
#' @export

hotspot_splot<-function(y1,y2,x,k,r,n,sz,l){
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
  #random generate number from negative binomial distribution where the mean parameter follows the pattern above
  generate<-function(s){
    g<-mean(rnbinom(30,mu=s,size=sz))
  }
  cresult<-sapply(result,generate)
  dresult<-matrix(cresult,nrow=2*x-1,ncol=2*x-1,byrow=T)

  #re-locate the pattern part on background and calibrate it to the index location as defined
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

  #draw the spot pattern
  graph<-filled.contour(x = 1:nrow(aresult),y = 1:ncol(aresult),
                        z = aresult, color.palette = myPalette,
                        plot.title = title(main = "Hotspot Spot-level Pattern ",
                                           xlab = "x-coordinate",ylab = "y-coordinate"),
                        plot.axes = {axis(1, seq(1, ncol(aresult), by = 5))
                          axis(2, seq(1, nrow(aresult), by = 5))},
                        key.title = title(main="Gene\n(counts)"),
                        key.axes = axis(4, seq(min(aresult), max(aresult), by = 20)),
  )
  finalresult<-list(aresult,graph)
  return(finalresult)
}
