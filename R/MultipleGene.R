install.packages("matrixcalc")
install.packages("lqmm")
library(lqmm)
library(matrixcalc)
library(MASS)
library(scatterpie)
library(RColorBrewer)
#not real-data based scRNA-seq data simulation

#N:cell number of each cell type
#k:number of different genes within each cell type 
#mean: kx1 verctor,mean of gene counts of each cell type
#cov: kxk matrix, covariance matrix
#disper: kx1 verctor of theta parameters in rnegbin function, variance=mu+mu^2/theta
#pref: name if the cell type

#name of cell type need to be defined outside the function
#if want to simulate real scRNA-seq data, dataset need to be set at first
prefix<-"Cell type A"
kernalmul1<-function(N,k,mean,cov,disper,status){
  if (status==0){
gene<-matrix(mvrnorm(n=1,mean,cov),nrow=1)
ge <- vector(mode = "list", length = k)
  for (i in 1:k){
    ge[[i]]<-matrix(rnbinom(N,mu=gene[i],size =disper[i]),ncol=1)
  }
suffix<-seq(1:N)
typelabel<-matrix(paste(prefix,suffix,sep="."),ncol=1)
type<-cbind(typelabel,ge[[1]])
for (j in 2:k){
  type<-cbind(type,ge[[j]])
 }
typef<-t(type)
}
  if (status==1){
    typef<-splatSimulate(splatEstimate(realdata))
  }
  return(typef)
}


#define distribution of locations and density
set.seed(42)
kernallode<-function(d){
uniform<-round(rgamma(n=1, shape=(d*0.8)/((d/0.3)/(d*0.8)), scale=(d/0.3)/(d*0.8)))
sparse<-round(rgamma(n=1, shape=(d*0.1)/((d/0.3)/(d*0.1)), scale=(d/0.3)/(d*0.1)))

#define distribution of density
highden1<-round(rgamma(n=uniform, shape=2.5, scale=1))
highden2<-round(rgamma(n=sparse, shape=2.5, scale=1))

lowden1<-round(rgamma(n=uniform, shape=0.8, scale=1))
lowden2<-round(rgamma(n=sparse, shape=0.8, scale=1))

#fix number location with cells out of all locations
loc1<-sort(sample(d,uniform,replace=F))
loc2<-sort(sample(d,sparse,replace=F))

#define cell number on locations
num1<-rpois(n=uniform,lambda=highden1) #uniform high
num2<-rpois(n=uniform,lambda=lowden1)  #uniform low
num3<-rpois(n=sparse,lambda=highden2)  #sparse high
num4<-rpois(n=sparse,lambda=lowden2)   #sparse low

  blank1<-matrix(c(rep(0,d)),ncol=1)
for (i in 1:d){
  blank1[loc1[i]]<-num1[i]
}

blank2<-matrix(c(rep(0,d)),ncol=1)
for (i in 1:d){
  blank2[loc1[i]]<-num2[i]
}

blank3<-matrix(c(rep(0,d)),ncol=1)
for (i in 1:d){
  blank3[loc2[i]]<-num3[i]
}

blank4<-matrix(c(rep(0,d)),ncol=1)
for (i in 1:d){
  blank4[loc2[i]]<-num4[i]
 }
return(list(blank1,blank2,blank3,blank4))
}

ld<-kernallode(400)

sum(ld[[1]]) #947 type A cell number in total
sum(ld[[2]]) #342 type B cell number in total
sum(ld[[3]]) #20 type C cell number in total
sum(ld[[4]]) #6 type D cell number in total
 
cellnum<-data.frame(xaxis=rep(1:20,20), yaxis=rep(1:20,each=20))
n<-nrow(cellnum)
cellnum$region <- factor(1:n)
cellnum$A <- as.numeric(ld[[1]])
cellnum$B <- as.numeric(ld[[2]])
cellnum$C <- as.numeric(ld[[3]])
cellnum$D <- as.numeric(ld[[4]])
cellnum$radius<- 0.31
#colnames(cellnum)[3+(1:length(num))]<-letters[1:length(num)]

p <- ggplot() + geom_scatterpie(aes(x=xaxis, y=yaxis, group=region, r=radius), data=cellnum,
                                cols=LETTERS[1:4], color=NA) + coord_equal()
p + geom_scatterpie_legend(cellnum$radius, x=0, y=0)
