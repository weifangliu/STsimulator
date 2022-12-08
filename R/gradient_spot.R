#Simulation of spot-level pattern (linear gradient)


#In this simulation model, gene expression counts of each spot (sum of gene counts of all cells within each spot) is modeled by negative binomial distribution. It allows to specify desired spatial linear-gradient-like patterns at spot level. One gene is considered in this simulation model


#' @param y2 index column location
#' @param y \code{2y-1} is the length of pattern
#' @param r change rate, the outer layer's mean is \code{1/r} times of that of the adjacent inner layer(column)
#' @param n mean parameter of the negative binomial distribution of index column
#' @param sz size parameter in alternative parametrization of negative binomial distribution
#' @param l length of background

#'
#' @return A graph with spatial linear-gradient pattern of SVG at spot level

#' @useDynlib sost
#' @importFrom gplots RColorBewer
#' @export

gradient_spot<-function(y2,y,r,k,n,sz,l){
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
    g<-mean(rnbinom(20,mu=s,size=sz))
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
  finalresult<-list(aresult,graph,dresult)
  return(finalresult)
}
