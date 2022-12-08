#Simulation of spot-level pattern (Streak)


#In this simulation model, gene expression counts of each spot (sum of gene counts of all cells within each spot) is modeled by negative binomial distribution. It allows to specify desired spatial streak-like patterns at spot level. One gene is considered in this simulation model


#' @param y1 one bound of index band
#' @param y2 another bound of index band

#'
#' @return A graph with spatial streak pattern of SVG at spot level

#' @useDynlib sost
#' @importFrom gplots RColorBewer
#' @export

streak_spot<-function(y1,y2,n,sz,l,base){
  result<-matrix(rep(0,l^2),nrow=l,ncol=l)
  result[,y1:y2]<-n
  result[,0:y1]<-base
  result[,y2:l]<-base
  generate<-function(s){
    g<-rnbinom(1,mu=s,size=sz)
  }
  cresult<-sapply(result,generate)
  aresult<-matrix(cresult,nrow=l,ncol=l,byrow=F)
  graph<-filled.contour(x = 1:ncol(aresult),y = 1:nrow(aresult),
                        z = t(aresult), color.palette = myPalette,
                        plot.title = title(main = "Streak Spot-level Pattern ",
                                           xlab = "x-coordinate",ylab = "y-coordinate"),
                        plot.axes = {axis(1, seq(1, ncol(aresult), by = 5))
                          axis(2, seq(1, nrow(aresult), by = 5))},
                        key.title = title(main="Gene\n(counts)"),
                        key.axes = axis(4, seq(min(aresult), max(aresult), by = 10)),
  )
  finalresult<-list(aresult,graph)
  return(finalresult)
}
