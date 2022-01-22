library("ESCO")
library("ellipsis")
setwd("/Users/ZhentaoYu/Desktop/ST/Simulation/mouse_ALM")
library(data.table)
library(dplyr)
library(Seurat)
library(patchwork)
geneinfo<-fread("mouse_ALM_2018-06-14_exon-matrix.csv",data.table = F)
#Normalization of count data 
geneinfo1<-NormalizeData(geneinfo, normalization.method = "LogNormalize", scale.factor = 10000)
cellinfo<-fread("mouse_ALM_2018-06-14_samples-columns.csv",data.table = F)
head(cellinfo)
summary(cellinfo$cluster)
#find top 2000 variable genes across all cells
geneinfo2<- FindVariableFeatures(geneinfo1, selection.method = "vst")
top<-geneinfo[order(geneinfo2$vst.variance.standardized,decreasing=T)[1:2000],]
Sncg<-cellinfo[cellinfo$subclass=="Sncg",]
Sst<-cellinfo[cellinfo$subclass=="Sst",]
Lamp5<-cellinfo[cellinfo$subclass=="Lamp5",]
Vip<-cellinfo[cellinfo$subclass=="Vip",]
Sncg1<-data.matrix(top[,colnames(top) %in% Sncg$sample_name])
Sst1<-data.matrix(top[,colnames(top) %in% Sst$sample_name])
Lamp51<-data.matrix(top[,colnames(top) %in% Lamp5$sample_name])
Vip1<-data.matrix(top[,colnames(top) %in% Vip$sample_name])
library(ESCO)
paramsSncg<-escoEstimate(Sncg1)
paramsSst<-escoEstimate(Sst1)
paramsLamp5<-escoEstimate(Lamp51)
paramesVip1<-escoEstimate(Vip1)
Sncgsim<-escoSimulateSingle(params=paramsSncg, verbose = T)
Sstsim<-escoSimulateSingle(params=paramsSst, verbose = T, nCells=2000)
Sstsim1<-escoSimulateSingle(params=paramsSst, verbose = T)
library(SingleCellExperiment)
Sncgcount<-counts(Sncgsim)
Sstcount<-counts(Sstsim)
#check zero proportion
sum(Sstcount==0)/sum(Sstcount!=0) #simulated data zero:non-zero: 133.7754
sum(Sst1==0)/sum(Sst1!=0) #oringinal data zero:non-zero: 18.121
#adjust sequencing depth before simulation: Sst1 example
nrow(Sst1) #2000. Sst1: gene by cell
#test on only changing sequecing depth
Sst2<-Sst1*5
Sstcount2<-counts(escoSimulateSingle(params=escoEstimate(Sst2), verbose = T))
sum(Sstcount2==0)/sum(Sstcount2!=0) #277.1101 even worse
#test on genes pick at prior
top2<-geneinfo[order(geneinfo2$vst.variance.standardized,decreasing=T)[c(1:100,30000:31000)],]
Sst3<-data.matrix(top2[,colnames(top2) %in% Sst$sample_name])
Sstcount3<-counts(escoSimulateSingle(params=escoEstimate(Sst3), verbose = T,nCells=6000))
sum(Sst3==0)/sum(Sst3!=0) #top 100 and 30000:31000 HVG zero:non-zero:1.78;  zero proportion:0.64
sum(Sstcount3==0)/sum(Sstcount3!=0)#simulated data 1.23; zero proportion: 0.55