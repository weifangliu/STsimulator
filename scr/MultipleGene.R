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

#mean function jiawen A type
mean<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  mean[i]<-round(runif(1,min=1000,max=1500))
}

#covariance functin jiawen A type
cov<-matrix(nrow=500,ncol=500)
diag(cov)<-runif(500,10,20)
for (i in 1:500){
for (j in min((i+1),500):500){
  cov[i,j]<-runif(1,min=-1,max=1)*sqrt(cov[i,i])*sqrt(cov[j,j])
  cov[j,i]<-cov[i,j]
 }
}
cov<-make.positive.definite(cov,tol=1*10^(-3))
#dispersion function jiawen A ype
disper<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  disper[i]<-round(runif(1,min=2,max=10))
}

jwa<-kernalmul1(1000,500,mean,cov,disper,0)
jwa[1:10,1:10]

#type B
prefix<-"Cell type B"
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

#mean function jiawen B type
mean<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  mean[i]<-round(runif(1,min=500,max=800))
}

#covariance functin jiawen B type
cov<-matrix(nrow=500,ncol=500)
diag(cov)<-runif(500,15,45)
for (i in 1:500){
  for (j in min((i+1),500):500){
    cov[i,j]<-runif(1,min=-1,max=1)*sqrt(cov[i,i])*sqrt(cov[j,j])
    cov[j,i]<-cov[i,j]
  }
}
cov<-make.positive.definite(cov,tol=1*10^(-3))
#dispersion function jiawen B ype
disper<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  disper[i]<-round(runif(1,min=20,max=30))
}

jwb<-kernalmul1(1000,500,mean,cov,disper,0)
jwb[1:10,1:10]

#type C
prefix<-"Cell type C"
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

#mean function jiawen C type
mean<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  mean[i]<-round(runif(1,min=3000,max=3200))
}

#covariance functin jiawen C type
cov<-matrix(nrow=500,ncol=500)
diag(cov)<-runif(500,5,6)
for (i in 1:500){
  for (j in min((i+1),500):500){
    cov[i,j]<-runif(1,min=-1,max=1)*sqrt(cov[i,i])*sqrt(cov[j,j])
    cov[j,i]<-cov[i,j]
  }
}
cov<-make.positive.definite(cov,tol=1*10^(-3))
#dispersion function jiawen C ype
disper<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  disper[i]<-round(runif(1,min=9,max=30))
}

jwc<-kernalmul1(1000,500,mean,cov,disper,0)
jwc[1:10,1:10]

#type D
prefix<-"Cell type D"
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

#mean function jiawen D type
mean<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  mean[i]<-round(runif(1,min=700,max=2000))
}

#covariance functin jiawen D type
cov<-matrix(nrow=500,ncol=500)
diag(cov)<-runif(500,40,60)
for (i in 1:500){
  for (j in min((i+1),500):500){
    cov[i,j]<-runif(1,min=-1,max=1)*sqrt(cov[i,i])*sqrt(cov[j,j])
    cov[j,i]<-cov[i,j]
  }
}
cov<-make.positive.definite(cov,tol=1*10^(-3))
#dispersion function jiawen D ype
disper<-matrix(rep(0,500),ncol=1)
for (i in 1:500){
  disper[i]<-round(runif(1,min=4,max=20))
}

jwd<-kernalmul1(1000,500,mean,cov,disper,0)
jwd[1:10,1:10]
#add rownames and colnames
# colnames(jwa)<-jwa[1,]
# jwa<-jwa[-1,]
# rownames(jwa)<-paste0("gene",1:500)
# 
# colnames(jwb)<-jwb[1,]
# jwb<-jwb[-1,]
# rownames(jwb)<-paste0("gene",1:500)
# 
# colnames(jwc)<-jwc[1,]
# jwc<-jwc[-1,]
# rownames(jwc)<-paste0("gene",1:500)
# 
# colnames(jwd)<-jwd[1,]
# jwd<-jwd[-1,]
# rownames(jwd)<-paste0("gene",1:500)

#four cell types, each has 500 cells
jwa1<-matrix(as.numeric(jwa),nrow=500,ncol=1000,byrow=F)
colnames(jwa1)<-paste0("Cell Type A.",1:1000)
rownames(jwa1)<-paste0("Gene",1:500)
jwb1<-matrix(as.numeric(jwb),nrow=500,ncol=1000,byrow=F)
colnames(jwb1)<-paste0("Cell Type B.",1:1000)
rownames(jwb1)<-paste0("Gene",1:500)
jwc1<-matrix(as.numeric(jwc),nrow=500,ncol=1000,byrow=F)
colnames(jwc1)<-paste0("Cell Type C.",1:1000)
rownames(jwc1)<-paste0("Gene",1:500)
jwd1<-matrix(as.numeric(jwd),nrow=500,ncol=1000,byrow=F)
colnames(jwd1)<-paste0("Cell Type D.",1:1000)
rownames(jwd1)<-paste0("Gene",1:500)


# test1<-kernalmul1(10,4,matrix(c(10,15,11,12),ncol=1),matrix(c(    1,  0.2, 0.2, 0.1,
#                                                                  0.2,  1, 0.5, 0.1,
#                                                                  0.2,0.5,   1, 0.3,
#                                                                 0.1,0.1, 0.3,   1), byrow=T, ncol=4),matrix(c(10,20,30,40),ncol=1),0)
# # # generate single cell data of cell type1
#
# type1_mean<-matrix(c(10,15,11,12),ncol=1)
# type1_gene<-matrix(mvrnorm(n=1,type1_logmean,type1_var),nrow=1)
# 
# type1_gene1<-matrix(rnegbin(5000,mu=type1_gene[1],theta=5),ncol=1)
# type1_gene2<-matrix(rnegbin(5000,mu=type1_gene[2],theta=4),ncol=1)
# type1_gene3<-matrix(rnegbin(5000,mu=type1_gene[3],theta=3),ncol=1)
# type1_gene4<-matrix(rnegbin(5000,mu=type1_gene[4],theta=2),ncol=1)
# prefix<-"type1"
# suffix<-seq(1:length(type1_gene1))
# type1_label<-matrix(paste(prefix,suffix,sep="."),ncol=1)
# type1<-cbind(type1_label,type1_gene1,type1_gene2,type1_gene3,type1_gene4)

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

#cells belonging to type A/B/C/D on the map
samplea<-sample(colnames(jwa1),sum(ld[[1]]))
sampleb<-sample(colnames(jwb1),sum(ld[[2]]))
samplec<-sample(colnames(jwc1),sum(ld[[3]]))
sampled<-sample(colnames(jwd1),sum(ld[[4]]))


lk<-ld

for (j in 1:4){
for (i in 1:400){
  lk[[j]][i]<-sum(ld[[j]][1:i])
 }
}

st<-matrix(nrow=500,ncol=400)
a<-matrix(nrow=500,ncol=400)
b<-matrix(nrow=500,ncol=400)
c<-matrix(nrow=500,ncol=400)
d<-matrix(nrow=500,ncol=400)
for (j in 1:500){
  st[j,1]<-sum(jwa1[j,which(colnames(jwa1)%in%samplea[1:(lk[[1]][1])])])+
    sum(jwb1[j,which(colnames(jwb1)%in%sampleb[1:(lk[[2]][1])])])+
    sum(jwc1[j,which(colnames(jwc1)%in%samplec[1:(lk[[3]][1])])])+
    sum(jwd1[j,which(colnames(jwd1)%in%sampled[1:(lk[[4]][1])])])
  for (i in 2:400){
  if (lk[[1]][i-1]<lk[[1]][i]){
      a[j,i]<-sum(jwa1[j,which(colnames(jwa1)%in%samplea[((lk[[1]][i-1])+1):(lk[[1]][i])])])
    } else{
      a[j,i]<-0
    }
    
  if (lk[[2]][i-1]<lk[[2]][i]){
      b[j,i]<-sum(jwb1[j,which(colnames(jwb1)%in%sampleb[((lk[[2]][i-1])+1):(lk[[2]][i])])])
    } else{
      b[j,i]<-0
    }
  if (lk[[3]][i-1]<lk[[3]][i]){
      c[j,i]<-sum(jwc1[j,which(colnames(jwc1)%in%samplec[((lk[[3]][i-1])+1):(lk[[3]][i])])])
    } else{
      c[j,i]<-0
    }
  if (lk[[4]][i-1]<lk[[4]][i]){
      d[j,i]<-sum(jwd1[j,which(colnames(jwd1)%in%sampled[((lk[[4]][i-1])+1):(lk[[4]][i])])])
    } else{
      d[j,i]<-0
    }
    st[j,i]<-a[j,i]+b[j,i]+c[j,i]+d[j,i]
  }
}
rownames(st)<-paste0("Gene.",1:500)
colnames(st)<-paste0("Spot.",1:400)

write.csv(st,"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulationST.csv")
write.csv(jwa1,"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_A_ref.csv")
write.csv(jwb1,"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_B_ref.csv")
write.csv(jwc1,"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_C_ref.csv")
write.csv(jwd1,"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_D_ref.csv")
write.csv(ld[[1]],"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_A_truth.csv")
write.csv(ld[[2]],"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_B_truth.csv")
write.csv(ld[[3]],"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_C_truth.csv")
write.csv(ld[[4]],"/Users/ZhentaoYu/Desktop/Bios992/JiawenSimulation_D_truth.csv")



#example
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


