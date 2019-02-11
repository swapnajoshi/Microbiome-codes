#setwd("C:/Users/swapnajoshi/Documents/UCLA_research/IBSmicrobiomeProject/biopsy_scor_microbiome/")
setwd("C:/Users/swapnajoshi-admin/Documents/swapna/")
install.packages("vegan",  dependencies=T)
install.packages("C:/Users/swapnajoshi-admin/Downloads/GO.db_2.0.2.zip", repos=NULL)
library(vegan)

# read OTU table

header1<-scan("otu_table_even_9188.biom.txt", nlines = 2, what = character())[8:91]
genus=read.table("otu_table_even_9188.biom.txt", sep="\t", skip=2, header=FALSE, row.names = 1)	
names(genus)<-header1

taxa.names = genus$ConsensusLineage

genus1 = as.matrix(genus[,-dim(genus)[2]])
dim(genus1)

# genus1 = scale(genus1, center=F, scale=colSums(genus1))
#genus1 = as.data.frame(t(genus1))
#colnames(genus1)<-taxa.names
save(genus1, file="R_output/genus1.rda")
save(genus, file="R_output/genus.rda")

#extract.name.level = function(x, level){
#  a=c(unlist(strsplit(x,';')),'Other')
#  paste(a[1:min(level,length(a))],collapse=';')
#}

#otu2taxonomy = function(x, level, taxa=NULL){
#  if(is.null(taxa))
#  {taxa = colnames(x)}
#  if(length(taxa)!=dim(x)[2])
#  {print("ERROR: taxonomy should have the same length
#         as the number of columns in OTU table")
#    return;}
#  level.names = sapply(as.character(taxa),
#                       function(x)
#                         extract.name.level(x,level=level))
#  t(apply(x, 1, function(y) tapply(y,level.names,sum)))
#}
 
# d.phylum = otu2taxonomy(genus1,level=8,taxa=taxa.names)


### mapping file
map1<-read.csv("ibs.mapping_new.csv",  header=T,row.names= 1, sep=',')
map2<-map1[row.names(map1)%in%colnames(genus1),] ## 2 samples did not have data

map2<-map2[colnames(genus1),]
match(row.names(map2),colnames(genus1))
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
# [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
# [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
# [76] 76 77 78 79 80 81 82 83
dim(genus1)
# [1] 83 4325
genus2<-t(genus1)
summary(colMeans(genus2))
hist(colMeans(genus2))
sorted.means<-rev(sort(colMeans((genus1))))
fix(sorted.means)

count.dat<-matrix(NA, nrow=dim(genus1)[1],ncol=dim(genus1)[2])
for (i in 1:dim(genus1)[1])
  for (j in 1: dim(genus1)[2])
  {
    if (genus1[i,j]>2)#> log(20,2) [1] 4.32
      
    { count.dat[i,j]<-1}
    
    else if (genus1[i,j]<=2)
    { count.dat[i,j]<-0}
  }

row.names(count.dat)<-row.names(genus1)
colnames(count.dat)<-colnames(genus1)

count.dat<-as.data.frame(count.dat)

count.dat$sum1<-apply(count.dat,1,sum)

##### Alteast 25% of the subjetcs have 25 counts or more per species-- 20% of 83 =16 samples
probesInt<-subset(count.dat, count.dat$sum1>=16)
ma.phylum<-genus1[row.names(genus1)%in%row.names(probesInt),]
dim(ma.phylum)
# 229 83


genus3<-log(ma.phylum+0.000001,2)
genus4<-apply(genus3,1, scale)
row.names(genus4)<-colnames(genus3)


genus4[1:5,1:5]
dim(genus4)
hc1<-hclust(dist(genus4), method="complete")
plot(hc1)

map<-read.delim("C:/Users/swapnajoshi-admin/Documents/Microbiome_IBS/ibs.mapping_new.csv", sep=",",row.names=1)
match(row.names(genus4),row.names(map))
map<-map[row.names(genus4),]
match(row.names(genus4),row.names(map))
labs=paste(row.names(genus4),map$Dx, sep="_")
row.names(genus4)<-labs
hc1<-hclust(dist(genus4), method="complete")
plot(hc1)

#otu6<-genus5[hc1$order,]
#colStates<-factor(row.names(otu6),levels=hc1$labels[hc1$order])
#Clust1<-cutree(hc1, k=7)
#table(Clust1)
#Clust2<-Clust1[hc1$order]

#only with Dx labels

labs=map$Dx
row.names(genus4)<-labs
hc1<-hclust(dist(genus4), method="complete")
plot(hc1)


Clust1<-cutree(hc1, k=5)
table(Clust1)
 1  2 3  4  5 
22  1 46 13  1 

otu6<-genus5[hc1$order,]
colStates<-factor(row.names(otu6),levels=hc1$labels[hc1$order])
Clust2<-Clust1[hc1$order]
#colClust<-factor(Clust2,levels=hc1$labels[hc1$order])
write.table(Clust2, file="ibs.clu.memb.heirrar_ave.csv", sep=",")
library(ggdendro)
library(ggplot2)
df2<-data.frame(Clust2)
head(df2)
rainbow(n=3)
[1] "#FF0000FF" "#00FF00FF" "#0000FFFF"

df2$color1<-NA
for ( i in 1: 52)
{if (df2[i,1]=="1")
	{df2[i,2]<-"#FF0000FF"}
else if (df2[i,1]=="2")
	{df2[i,2]<-"#00FF00FF"}
else if (df2[i,1]=="3")
	{df2[i,2]<-"#0000FFFF"}
else if (df2[i,1]=="4")
	{df2[i,2]<-"#0000FFFF"}
}
clust2<-as.data.frame(Clust2)
df.hc<-cbind(row.names(map[map$Dx=="HC",]), c(rep("HC",31)))
row.names(df.hc)<-df.hc[,1]
df.hc<-as.data.frame(df.hc[,-1])
colnames(df.hc)<-colnames(clust2)
clust3<-rbind(clust2,df.hc)

map1<-merge(map, clust3, by="row.names")
row.names(map1)<-map1[,1]
map1$cluster1<-NULL
map1$Row.names<-NULL

#write the membership to mapping file


map1$cluster1<-NA
for( i in 1: dim(clust3)[1])
{
if(map1[i,79]==1){map1[i,80]<-1}
else if(map1[i,79]!=1) {map1[i,80]<-0}
}

map1$cluster2<-NA
for( i in 1: dim(df2)[1])
{
if(map1[i,79]==2){map1[i,81]<-1}
else if(map1[i,79]!=2) {map1[i,81]<-0}
}

map1$cluster3<-NA
for( i in 1: dim(df2)[1])
{
if(map1[i,79]==3 |map1[i,79]==4){map1[i,82]<-1}
else if(map1[i,79]!=3) {map1[i,82]<-0}
}

map1$HC<-NA
for( i in 1: dim(df2)[1])
{
if(map1[i,79]=="HC"){map1[i,83]<-1}
else if(map1[i,79]!="HC") {map1[i,83]<-0}
}
## include member 4 in 3 
map1$Clust2<-gsub("4","3",map1$Clust2)
map2<-map1[,c(1:77,79:83,78)]
colnames(map2)[78]<-"cluster"


write.table(map2, file="bigIbsMapping.csv",sep=",", col.names=NA)
colnames(map2)
littleIbsMapping<-map2[,c(1,2,3,4,5,6,8,9,10,30:42,76:83)]
write.table(littleIbsMapping, file="littleIbsMapping.csv",sep=",", col.names=NA)

#d.phylum<-t(d.phylum)

#count.dat<-matrix(NA, nrow=dim(d.phylum)[1],ncol=dim(d.phylum)[2])
#for (i in 1:dim(d.phylum)[1])
# for (j in 1: dim(d.phylum)[2])
#  {
#    if (d.phylum[i,j]>25)#> log(20,2) [1] 4.32
#      
#    { count.dat[i,j]<-1}
#    
#    else if (d.phylum[i,j]<=25)
#    { count.dat[i,j]<-0}
#  }
#
#row.names(count.dat)<-row.names(d.phylum)
#colnames(count.dat)<-colnames(d.phylum)
#
#count.dat<-as.data.frame(count.dat)
#
#count.dat$sum1<-apply(count.dat,1,sum)
#

##### Alteast 25% of the subjetcs have 25 counts or more per species-- 25% of 83 =20 samples
#probesInt<-subset(count.dat, count.dat$sum1>=20)
#ma.phylum<-d.phylum[row.names(d.phylum)%in%row.names(probesInt),]


ma.phylum<-genus1[,names(sorted.means[which(sorted.means>=10)])]
dim(ma.phylum)
# 83 71
vars<-row.names(ma.phylum)
var.dat<-ma.phylum[row.names(map2),]
match(rownames(var.dat), row.names(map2))
par(mar = rep(2, 4))
boxplot(var.dat[,1], map2$Dx)

var.dat.spec<-cbind(var.dat,map2)
save(genus1, ma.phylum, map2, vars, var.dat, genus, file="R_output/OtuToDataFrame.Rda")

save(ma.phylum, file="R_output/sel.bact.phy.matrix.rda")

# Alpha diversitites

num_genera=specnumber(ma.phylum)
chao1=estimateR(ma.phylum*100000000)[2,]
shannon=diversity(ma.phylum,  "shannon")

save(shannon, file="shannon.rda")

# plot alpha diversity

par(mar = rep(2, 4))
boxplot(num_genera  ~	var.dat.spec[,78], col=rainbow(3), main="Number of genera")
boxplot(chao1~  var.dat.spec[,78],  col=rainbow(3),	 main="Chao1  index")
boxplot(shannon~  var.dat.spec[,78],	  col=rainbow(3),  main="Shannon  evenness")


 summary(aov(num_genera  ~  var.dat.spec[,96]))
#                   Df Sum Sq Mean Sq F value Pr(>F)  
#var.dat.spec[, 96]  3    801  267.17    3.56 0.0179 *
#Residuals          79   5929   75.05    

# Comparison of bacterial communities
# 1. which species are different between ibs and hc,
# 2. diffrence between bowel habbit subtypes

source("https://bioconductor.org/biocLite.R")
biocLite("source("https://bioconductor.org/biocLite.R")
biocLite("dplyr")
# heatmap

#########################################################
#setwd("C:/Users/swapnajoshi/Documents/UCLA_research/IBSmicrobiomeProject/biopsy_scor_microbiome/")
load("OtuToDataFrame.Rda")
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
#WGCNA
library(WGCNA)
library(flashClust)
#Heatmap
library(dplyr)
library(NMF)
library(RColorBrewer)

# trait data 
traitData = map2 #from above

dim(ma.phylum)
ma.phylum[1:2,1:2]
#fix(d.phylum)
# datExpr0<-ma.phylum
#datExpr0<-subset(ma.phylum, map2$Dx=="IBS")
#map2.ibs<-subset(map2, map2$Dx=="IBS")
# row.names(datExpr0)<-substr(row.names(datExpr0),8,13)


#gsg <- goodSamplesGenes(datExpr0, verbose = 3)
#gsg$allOK




#sampleTree = flashClust(dist(datExpr0), method = "average")
#match(row.names(map2.ibs),row.names(datExpr0)) ## ok
#labels1=map2.ibs$Dx
#dat.ordered<-datExpr0[sampleTree$order,]

# otuData.ibs<-subset(datExpr0, row.names(datExpr0)%in% [traitData grep("IBS",colnames(otuData)),])
# otuData.hc<-datExpr0[grep("HC",colnames(otuData)),]

#sizeGrWindow(12,9)
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#plot(sampleTree, main = "Sample clustering to detect outliers_ibs_new", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
#abline(h = 100000, col = "red")
#load("memb1.rda")
############################################# hc and ibs together 
#clust1 = cutreeStatic(sampleTree, cutHeight = 100000, minSize = 1)
#table(clust1)
# clust1
# 1  2  3  4   
# 25 18  6  3   

#memb1<-as.matrix(cbind(rownames(datExpr0),clust1))
#memb1<-as.data.frame(memb1)
#memb1$labels<-labels1

#   HC IBS
# 1 15  25
# 2 12  18
# 3  2   6
# 4  1   3
# 5  1   0

# memb1.ibs<-subset(memb1, memb1$labels=="IBS")
# memb1.hc<-subset(memb1, memb1$labels=="HC")
# table(memb1$clust1,memb1$labels)
#save(memb1,file="memb1.rda")
# datExpr1<-t(apply(datExpr0,1,scale))
# row.names(datExpr1)<-row.names(datExpr0)
# colnames(datExpr1)<-colnames(datExpr0)


exp.clu1<-t(ma.phylum[row.names(ma.phylum)%in%row.names(map2[map2[,25]=="1",]),])
exp.clu2<-t(ma.phylum[row.names(ma.phylum)%in%row.names(map2[map2[,25]=="2",]),])
exp.clu3<-t(ma.phylum[row.names(ma.phylum)%in%row.names(map2[map2[,25]=="3",]),])
exp.hc<-t(ma.phylum[row.names(ma.phylum)%in%row.names(map2[map2[,25]=="HC",]),])
# 
save(exp.clu1, file="R_output/exp.clu1.rda")
save(exp.clu2, file="R_output/exp.clu2.rda")
save(exp.clu3, file="R_output/exp.clu3.rda")
save(exp.hc, file="R_output/exp.hc.rda")

########################################################################

# hc<-t(datExpr1[row.names(datExpr1)%in%memb1.hc[,1],])
#ibs.pre<-exp.clu1
# ibs.pre1<-ibs.pre[,!colnames(ibs.pre)%in%colnames(hc)]
# dim(ibs.pre1)
# # [1]  8 36
#ibs.norLike<-exp.clu2
# dim(ibs.norLike)
# # [1]  8 16
# dim(hc)
# [1]  8 31
# otuData.ibs<-otuData.zscore[,grep("IBS",colnames(otuData.zscore))]
# otuData.hc<-otuData.zscore[,grep("HC",colnames(otuData.zscore))]
# sex<-function(mol.biol) {
#   if (mol.biol=="1") "blue" #red
#   else if (mol.biol=="2") "pink" #yellow
#   else if (mol.biol=="NA") "grey" #Black
# }
# 
# BowelHabit<-function(mol.biol) {
#   if (mol.biol=="2") "blue" #red
#   else if (mol.biol=="1") "green" #yellow
#   else if (mol.biol=="4") "red" #Black
#   else if (mol.biol=="6") "yellow" #Black
#   else if (mol.biol=="5") "yellow" #Black
#   else if (mol.biol=="NA") "grey"
# }
# map3<-map2
# map3$id<-row.names(map3)
# map3<-map3[,c(44,1:43)]
# cc.sex<-apply(as.matrix(colnames(ibs.pre1)[1:36]),1,vlookup,map3,10)
# cc1<-lapply(cc.sex, as.character)
# cc.col1<-unlist(lapply(cc.sex, sex))
# 
# cc.bh<-apply(as.matrix(colnames(ibs.pre1)[1:36]),1,vlookup,map3,9)
# cc.col2<-unlist(lapply(cc.bh, BowelHabit))
# col.mat<-cbind(cc.col1,cc.col2)
# 
# col.j<-jet.colors(75)
# 
# 
# 
# png("hm.ibs.pre.png", bg="white", res=300, width=3000, height=3000)
# hv2<-heatmap.3(as.matrix(ibs.pre1),na.rm=TRUE,scale="none",  
#                #		RowSideColor=col.mat,
#                ColSideColors=col.mat,
#                #		col=col.j,
#                #		zlim=c(0,1),
#                Colv=T,Rowv=T,
#                cexRow=1,cexCol=1,
#                #		hclustfun = dist1,
#                # keep.dendro = NA,
#                main = paste("IBS_HC", dim(ibs.pre1)[1],"Probes; ",dim(ibs.pre1)[2],"Samples"),
#                labRow=substr(row.names(ibs.pre1),17,30)
# )
# dev.off()
# row.order <- rev(hv2$rowInd)
# #		col.order <- colnames(tum.dic.pr1[ ,hv1$colInd])
# ibs.norLike<-ibs.norLike[row.order,]
# 
# cc.sex<-apply(as.matrix(colnames(ibs.norLike)[1:16]),1,vlookup,map3,10)
# cc1<-lapply(cc.sex, as.character)
# cc.col1<-unlist(lapply(cc.sex, sex))
# 
# cc.bh<-apply(as.matrix(colnames(ibs.norLike)[1:16]),1,vlookup,map3,9)
# cc.col2<-unlist(lapply(cc.bh, BowelHabit))
# col.mat<-cbind(cc.col1,cc.col2)
# 
# 
# png("hm.ibs.normalLike.png", bg="white", res=300, width=3000, height=3000)
# hv2<-heatmap.3(as.matrix(ibs.norLike),na.rm=TRUE,scale="none",  
#                #		RowSideColor=col.mat,
#                ColSideColors=col.mat,
#                #		col=col.j,
#                #		zlim=c(0,1),
#                Colv=T,Rowv=NA,
#                cexRow=1,cexCol=1,
#                #		hclustfun = dist1,
#                # keep.dendro = NA,
#                main = paste("IBS", dim(ibs.norLike)[1],"Probes; ",dim(ibs.norLike)[2],"Samples"),
#                labRow=substr(row.names(ibs.norLike),17,30)
# )
# dev.off()
# 
# # row.order <- rev(hv2$rowInd)
# #		col.order <- colnames(tum.dic.pr1[ ,hv1$colInd])
# hc<-hc[row.order,]
# 
# cc.sex<-apply(as.matrix(colnames(hc)[1:31]),1,vlookup,map3,10)
# cc1<-lapply(cc.sex, as.character)
# cc.col1<-unlist(lapply(cc.sex, sex))
# col.mat<-cbind(cc.col1,cc.col1)
# 
# png("hm.hc.png", bg="white", res=300, width=3000, height=3000)
# hv2<-heatmap.3(as.matrix(hc),na.rm=TRUE,scale="none",  
#                #		RowSideColor=col.mat,
#                ColSideColors=col.mat,
#                #		col=col.j,
#                #		zlim=c(0,1),
#                Colv=T,Rowv=NA,
#                cexRow=1,cexCol=1,
#                #		hclustfun = dist1,
#                # keep.dendro = NA,
#                main = paste("hc", dim(hc)[1],"Probes; ",dim(hc)[2],"Samples"),
#                labRow=substr(row.names(hc),4,15)
# )
# dev.off()

####################################################################################

# heatmap for reference
#scale together

z<-apply(ma.phylum,1,scale)
row.names(z)<-colnames(ma.phylum)
colnames(z)<-row.names(ma.phylum)
clu1.z<-z[,colnames(z)%in%colnames(exp.clu1)]
clu2.z<-z[,colnames(z)%in%colnames(exp.clu2)]
clu3.z<-z[,colnames(z)%in%colnames(exp.clu3)]
hc.z<-z[,colnames(z)%in%colnames(exp.hc)]

# dat.ordered.ibs<-dat.ordered[row.names(dat.ordered)%in%row.names(map2.ibs[map2.ibs$Dx=="IBS",]),]

# datExpr1<-apply(dat.ordered.ibs,1,scale)
# row.names(datExpr1)<-colnames(dat.ordered.ibs)
# colnames(datExpr1)<-row.names(dat.ordered.ibs)
# dist1=hclust(dist(datExpr1), method="average")
library(heatmap.plus)
library(matlab)
col.j<-jet.colors(75)
png("R_output/hm.clu1.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(clu1.z),na.rm=TRUE,scale="none",  
               #		RowSideColor=col.mat,
               # ColSideColors=col.mat,
               		col=col.j,
               #		zlim=c(0,1),
               Colv=T,Rowv=T,
               # cexRow=1,cexCol=1,
               # distfun = dist1,
               keep.dendro = TRUE,
               main = paste("IBS_Predominant", dim(clu1.z)[1],"Probes; ",dim(clu1.z)[2],"Samples"),
               labRow=row.names(clu1.z))

dev.off()
row.order <- rev(hv2$rowInd)
ibs.norLike.z<-ibs.norLike.z[row.order,]


png("hm.ibs.normal-like.all.taxa_new.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(ibs.norLike.z),na.rm=TRUE,scale="none",  
                  #		RowSideColor=col.mat,
                  # ColSideColors=col.mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=T,Rowv=NA,
                  # cexRow=1,cexCol=1,
                  # distfun = dist1,
                  keep.dendro = TRUE,
                  main = paste("Normal-like IBS", dim(ibs.norLike.z)[1],"Probes; ",dim(ibs.norLike.z)[2],"Samples"),
                  labRow=row.names(ibs.norLike.z))

dev.off()


# row.order <- rev(hv2$rowInd)
hc.z<-hc.z[row.order,]


png("hm.hc.alltaxa_new.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(hc.z),na.rm=TRUE,scale="none",  
                  #		RowSideColor=col.mat,
                  # ColSideColors=col.mat,
                  col=col.j,
                  #		zlim=c(0,1),
                  Colv=T,Rowv=NA,
                  # cexRow=1,cexCol=1,
                  # distfun = dist1,
                  keep.dendro = TRUE,
                  main = paste("Healthy controls", dim(hc.z)[1],"Probes; ",dim(hc.z)[2],"Samples"),
                  labRow=row.names(hc.z))

dev.off()


# pie chart

# ibs.pre,ibs.norLike,exp.hc


# slices <- c(10, 12, 4, 16, 8)
png("pie.chart.all.taxa_new.png", res=300, width=3000, height=1500)

par(mfrow=  c(1,3), mar=  c(1, 2, 2, 8))
slices.ibs.pre<-round(rowMeans(ibs.pre), digits=0)



pct <- round(slices.ibs.pre/sum(slices.ibs.pre)*100)
pct1<-subset(pct,pct>2)
lbls <- names(pct1)
# lbls <- paste(lbls, pct) # add percents to labels
# lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(pct1,labels = NA, col=rainbow(length(lbls)),radius = 1.0,
    main="Pie Chart IBS") 



slices.ibs.norLike<-subset(round(rowMeans(ibs.norLike), digits=0), names(round(rowMeans(ibs.norLike), digits=0))%in%lbls)
lbls
# slices <- c(10, 12, 4, 16, 8)


# lbls <- names(slices.ibs.norLike)
pct <- round(slices.ibs.norLike/sum(slices.ibs.pre)*100)
# lbls <- paste(lbls, pct) # add percents to labels
# lbls <- paste(lbls,"%",sep="") # ad % to labels
# pct1<-subset(pct,pct>2)
# lbls <- names(pct1)
pie(pct,labels = NA, col=rainbow(length(lbls)),radius = 1.0,
    main="Pie Chart Normal-like IBS") 



slices.hc<-subset(round(colMeans(exp.hc), digits=0), names(round(colMeans(exp.hc), digits=0))%in%lbls)

# slices.hc<-round(colMeans(exp.hc[,1:7]), digits=0)

# slices <- c(10, 12, 4, 16, 8)

  
  # lbls <- names(slices.hc)
pct <- round(slices.hc/sum(slices.hc)*100)
# lbls <- paste(lbls, pct) # add percents to labels
# lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(pct,labels = NA, col=rainbow(length(lbls)),radius = 1.0,
    main="Pie Chart healthy controls") 

legend("topright", c("k__Bacteria; p__Actinobacteria",  "k__Bacteria; p__Bacteroidetes" ,  "k__Bacteria; p__Cyanobacteria" , 
                    "k__Bacteria; p__Firmicutes"  ,    "k__Bacteria; p__Proteobacteria",  "k__Bacteria; p__Tenericutes" ,   
                     "k__Bacteria; p__Verrucomicrobia"), cex=1, fill=rainbow(length(lbls)))
dev.off()




clust = cutreeStatic(sampleTree, cutHeight = 120000, minSize = 80)
table(clust)
# clust 0 contains the samples we want to keep.
keepSamples = (clust==0)

datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData<-map2.ibs

dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(1,2,3,4,5,6,7,10, 12,13,14, 15,16,17,18,19,20,43)];
# allTraits = allTraits[, c(2, 8:56) ];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
ibsSamples = rownames(datExpr);
traitRows = match(ibsSamples, row.names(allTraits));
datTraits = allTraits[traitRows, ];
match(ibsSamples, row.names(datTraits));
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 

# rownames(datTraits) = allTraits[traitRows, ];
collectGarbage();

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
datTraits<-as.data.frame(datTraits)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
png("traitsSampleTree_ibs_new.png", width=6500,height=5500, res=600)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(datTraits,sampleTree2,traitColors,datExpr,ibsSamples,allTraits,datExpr0,sampleTree, file="wgcna.part1.files.RData")

########### Do the groups exist in HCs?
load("otuData.ibs.rda")
# exp.hc.sel<-exp.hc[row.names(exp.hc)%in%row.names(exp.ibs.sel),]
otuData = as.data.frame(t(otuData.hc))
dim(otuData)
# [1] 226  31
names(otuData)
fix(otuData)

datExpr0 = as.data.frame(t(otuData))
# names(datExpr0) = otuData$ProbeID
# rownames(datExpr0) = names(otuData)[-c(1:8)]
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# [1] TRUE
sampleTree = flashClust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
png("SampleClusteringOutlierDectection_hc.png", width=2300, height=1800, res=200)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
# Plot a line to show the cut
abline(h = 49, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 49, minSize = 50)
table(clust)
# clust
# 0 
# 15 
# clust 0 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData = read.csv("sjoshi_samples_miRNA_colon_020116.csv")
row.names(traitData)<-traitData[,2]
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(1,2,3,4,5,6,7,10,12:20,43,44:67)];
# allTraits = allTraits[, c(2, 8:56) ];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
ibsSamples = rownames(datExpr);
traitRows = match(ibsSamples, row.names(allTraits));
datTraits = allTraits[traitRows, ];
match(ibsSamples, row.names(datTraits));
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15

# rownames(datTraits) = allTraits[traitRows, ];
collectGarbage();

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
datTraits<-as.data.frame(datTraits)
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
png("traitsSampleTree_hc.png", width=6500,height=5500, res=600)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

######################## split patients into 2 groups

clust1 = cutreeStatic(sampleTree2, cutHeight = 29, minSize = 1)
table(clust1)

memb1<-as.matrix(cbind(row.names(datExpr),clust1))
exp.clu1<-t(datExpr[row.names(datExpr)%in%memb1[memb1[,2]=="1",1],])
exp.clu2<-t(datExpr[row.names(datExpr)%in%memb1[memb1[,2]=="2",1],])

save(memb1,file="memb1.rda")

save(exp.clu1, file="exp.clu1.rda")
save(exp.clu2, file="exp.clu2.rda")




# gene expression differnces between cluster 1 and cluster 2
# asociation with traits

# t-test for differnces between two groups







####################################### Constructing adjacency matrix
load("wgcna.part1.files.RData")


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes


powers=c(seq(1,10,by=1),seq(12,20,by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(10, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.2,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 5;
adjacency1 = adjacency(datExpr, power = softPower);

TOM = TOMsimilarity(adjacency1);
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");

# sizeGrWindow(12,9)
png("Gene clustering on TOM-based dissimilarity.png", width=2500, height=1800, res=300)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

minModuleSize = 1;

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicMods
# 2 2 3 2       3 
# 1 1 2 1 0 0 0 2 

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# dynamicColors
# blue      grey turquoise 
# 2         3         3 

sizeGrWindow(8,6)

png("gene dendrogram module colors.png",width=3000, height=1500, res=300 )

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# signif(cor(MEs,use="p"),2)

MEDiss = 1-cor(MEs);

METree = hclust(as.dist(MEDiss), method = "average");
# sizeGrWindow(7, 6)
png("Clustering of module eigengenes.png", width=3000, height=1500, res=300)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

MEDissThres = 0.25

save(datExpr, dynamicColors, METree, MEs, dynamicColors, MEDiss, geneTree, MEList, file="modules.part2.RData")
# abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors;

mergedMEs = merge$newMEs;

# sizeGrWindow(12, 9)
png("Dynamic Tree Cut_merged dynamic.png", width=3000, height=1500, res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors
save(MEs, moduleLabels, moduleColors, geneTree, file = "ibs_module_data_files.RData")

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);


MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
signif(cor(MEs, use="p"), 2)
> # MEturquoise MEbrown MEyellow MEblue MEgreen MEblack    MEred MEgrey
  > # MEturquoise       1.000  -0.670 -0.65000  -0.53   -0.79   -0.39 -0.39000  0.076
  > # MEbrown          -0.670   1.000  0.72000   0.16    0.63   -0.13  0.12000 -0.095
  > # MEyellow         -0.650   0.720  1.00000   0.64    0.68    0.40  0.00034  0.030
  > # MEblue           -0.530   0.160  0.64000   1.00    0.70    0.59 -0.29000  0.100
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.table(moduleTraitPvalue, file="moduleTraitPValue.csv", sep=",")
anySig<-apply(moduleTraitPvalue<0.05,2,any)




# sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
png("module-trait relationships.png",width=2800,height=2200, res=300 )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

hierTOM = hclust( as.dist(dissTOM),method="average");
#The reader should vary the height cut -  off parameter h1
#(related to the y-axis of dendrogram)in the following 
colorStaticTOM =  as.character(cutreeStaticColor(hierTOM, cutHeight=1, minSize=20))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                    deepSplit=2, pamRespectsDendro = FALSE))

sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,    colors  = data.frame(colorStaticTOM,colorDynamicTOM, colorDynamicHybridTOM),
                    dendroLabels = FALSE, marAll =  c(1, 8, 3, 1),  main ="Gene dendrogram and module colors, TOM dissimilarity")


datME=moduleEigengenes(datExpr,dynamicColors)$eigengenes
signif(cor(datME, use="p"), 2)
# MEblue MEgrey MEturquoise
# MEblue       1.000  0.015       0.210
# MEgrey       0.015  1.000      -0.088
# MEturquoise  0.210 -0.088       1.000

dissimME=(1-t (cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average")
# Plot the eigengene dendrogram par (mfrow= c(1,1))
png("Clustering tree based of the module eigengenes.png")
plot (hclustdatME, main="Clustering tree based of the module eigengenes")
dev.off()

png("MEpairs.bloating.png", res=200, width=1500, height=1000)
# sizeGrWindow(8,9)
datTraits<-as.data.frame(datTraits)
datTraits$BSQ_Bloating[which(is.na(datTraits$BSQ_Bloating))]<-0
y= datTraits$BSQ_Bloating
plotMEpairs(datME, y=y)
dev.off()


# heatmap plots of module expression

sizeGrWindow(8,9)
par(mfrow=  c(3,1), mar=  c(1, 2, 4, 1))
which.module=  "turquoise";
# png("turquise_heatmap.png", width=1800, height= 1500, res=300)
plotMat(t(scale (datExpr[,dynamicColors== which.module ]) ),nrgcols=30,rlabels=T, 
        clabels=T,rcols= which.module,title=  which.module )
# dev.off()

which.module=  "blue";
# png("turquise_heatmap.png", width=1800, height= 1500, res=300)
plotMat(t(scale (datExpr[,dynamicColors== which.module ]) ),nrgcols=30,rlabels=T, 
        clabels=T,rcols= which.module,title=  which.module )
# dev.off()

which.module=  "grey";
# png("turquise_heatmap.png", width=1800, height= 1500, res=300)
plotMat(t(scale (datExpr[,dynamicColors== which.module ]) ),nrgcols=30,rlabels=T, 
        clabels=T,rcols= which.module,title=  which.module )
# dev.off()



sizeGrWindow(8,7);
which.module="blue"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1) , mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")


# datTraits<-as.data.frame(datTraits)
# datTraits$BSQ_UsualSeverity[which(is.na(datTraits$BSQ_UsualSeverity))]<-0



#####################################################
y= datTraits$BSQ_UsualSeverity
signif(cor(y,datME, use="p"),2)
#       MEblack MEblue MEbrown MEgreen MEgrey MEred MEturquoise MEyellow
# [1,]    0.28  -0.02   0.054   0.012  -0.19  0.54       -0.31      0.1
## only red is significant

cor.test(y, datME$MEred)[[3]]
# 0.003215997


sizeGrWindow(8,7);
which.module="red"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1) , mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="Red module", cex.main=2,
        ylab="eigengene expression",xlab="array sample")


GS1<-as.numeric(cor(y,datExpr,use="p"))
GeneSignificance<-abs(GS1)
ModuleSignificance<-tapply(GeneSignificance,dynamicColors,mean,na.rm=T)

#   black       blue      brown      green       grey        red  turquoise     yellow 
# 0.21962815 0.11689723 0.11305903 0.09067409 0.18572230 0.38949867 0.22747041 0.11726103 

sizeGrWindow(8,7)
par(mfrow(c(1,1)))
png("moduleSignificance.png")
plotModuleSignificance(GeneSignificance, dynamicColors)
dev.off()



severity = as.data.frame(datTraits$BSQ_UsualSeverity);
names(severity) = "severity"

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


geneTraitSignificance = as.data.frame(cor(datExpr, severity, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(severity), sep="");
names(GSPvalue) = paste("p.GS.", names(severity), sep="");
sigGSPvalue<-subset(GSPvalue, GSPvalue$p.GS.severity<0.05)
# annotate the probes
# annot<-read.delim("HTA-2_0.na33.1.hg19.transcript_gene_level_annotations.csv",sep=",")
# row.names(annot)<-annot[,1]
# gender.annot<-annot[row.names(annot)%in%row.names(sigGSPvalue),]

# correlation plot between gender and genes 

severity.genes<-datExpr[,colnames(datExpr)%in%row.names(sigGSPvalue)]

severity1 = as.data.frame(cbind(row.names(datTraits),datTraits$BSQ_UsualSeverity))
match(row.names(severity.genes), severity1[,1])
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34

df1<-cbind(severity.genes,severity1)
df1$V1<-NULL
cor.test(as.numeric(df1[,1]),as.numeric(df1[,15]), use='pairwise.complete.obs')[[3]]

cor.p<-matrix(NA, nrow=15, ncol=2)
row.names(cor.p)<-colnames(df1)
colnames(cor.p)<-c("cor1","pvalue")

for (i in 1:15)
{
  cor.p[i,1]<-cor.test(as.numeric(df1[,i]),as.numeric(df1[,15]), use='pairwise.complete.obs')[[4]]
  cor.p[i,2]<-cor.test(as.numeric(df1[,i]),as.numeric(df1[,15]), use='pairwise.complete.obs')[[3]]
}

# for ( i in 1:156)

#   
# {plotM<-function(x=1:156, y=157) plot

for (i in 1:15)
{
  png(paste(i,".png", sep=""))
  plot(df1[,i], df1[,15])
  dev.off()
}

#box with jitter
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNA_new_analyses_IBS9232015/wgcna_miRNA")
load("wgcna.part1.files.RData")


colnames(df1)<-gsub("\\|0.007","",colnames(df1))
colnames(df1)<-gsub("\\|0.033","",colnames(df1))
colnames(df1)<-gsub("\\|0","",colnames(df1))
colnames(df1)<-gsub("-","_",colnames(df1))

colnames(df1)

# [1] "hsa_miR_302d_3p" "hsa_miR_631"     "hsa_miR_19a_3p"  "hsa_miR_146b_5p" "hsa_miR_4516"    "hsa_miR_19b_3p" 
# [7] "hsa_miR_378e"    "hsa_miR_1283"    "hsa_miR_324_5p"  "hsa_miR_296_5p"  "hsa_miR_98"      "hsa_miR_1253"   
# [13] "hsa_miR_612"     "hsa_miR_188_5p"  "V2"  

ggplot(na.omit(df1), aes(x=V2,y=hsa_miR_302d_3p))+geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", outlier.shape = 16, outlier.size = 2, notch = F, notchwidth = 0.5)
bar1<-qplot(factor(V2),hsa_miR_302d_3p,data=na.omit(df1), geom=c("boxplot","jitter"), fill=factor(V2), main="BSQ_severity vs hsa-miR-302d-3p",xlab= "BSQ_Severity", ylab="hsa-miR-302d-3p")
x<-bar1 + stat_boxplot(geom ='errorbar')
ggsave(bar1, file="hsa-miR-302d-3p.boxplot.png", width = 8, height = 5, units = "in", dpi = 200)


# remove na from 
df1<-df1[-9,]


for (i in 1:15)
{
  
  x<-qplot(factor(df1[,15]),df1[,i], geom=c("boxplot","jitter"), fill=factor(df1[,15]), main="BSQ_severity",xlab= "BSQ_Severity") 
  y<-x + stat_boxplot(geom ='errorbar')
  print(y)
  ggsave(paste(i,"boxplot.png", sep=""))
}




################################################################################


#####################################################
y= datTraits$BSQ_UsualSeverity
signif(cor(y,datME, use="p"),2)
#       MEblack MEblue MEbrown MEgreen MEgrey MEred MEturquoise MEyellow
# [1,]    0.28  -0.02   0.054   0.012  -0.19  0.54       -0.31      0.1
## only red is significant

cor.test(y, datME$MEred)[[3]]
# 0.003215997


sizeGrWindow(8,7);
which.module="red"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1) , mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,dynamicColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="Red module", cex.main=2,
        ylab="eigengene expression",xlab="array sample")


GS1<-as.numeric(cor(y,datExpr,use="p"))
GeneSignificance<-abs(GS1)
ModuleSignificance<-tapply(GeneSignificance,dynamicColors,mean,na.rm=T)

#   black       blue      brown      green       grey        red  turquoise     yellow 
# 0.21962815 0.11689723 0.11305903 0.09067409 0.18572230 0.38949867 0.22747041 0.11726103 

sizeGrWindow(8,7)
par(mfrow(c(1,1)))
png("moduleSignificance.png")
plotModuleSignificance(GeneSignificance, dynamicColors)
dev.off()



severity = as.data.frame(datTraits$BSQ_UsualSeverity);
names(severity) = "severity"

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


geneTraitSignificance = as.data.frame(cor(datExpr, severity, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(severity), sep="");
names(GSPvalue) = paste("p.GS.", names(severity), sep="");
sigGSPvalue<-subset(GSPvalue, GSPvalue$p.GS.severity<0.05)
# annotate the probes
# annot<-read.delim("HTA-2_0.na33.1.hg19.transcript_gene_level_annotations.csv",sep=",")
# row.names(annot)<-annot[,1]
# gender.annot<-annot[row.names(annot)%in%row.names(sigGSPvalue),]

# correlation plot between gender and genes 

severity.genes<-datExpr[,colnames(datExpr)%in%row.names(sigGSPvalue)]

severity1 = as.data.frame(cbind(row.names(datTraits),datTraits$BSQ_UsualSeverity))
match(row.names(severity.genes), severity1[,1])
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34

df1<-cbind(severity.genes,severity1)
df1$V1<-NULL
cor.test(as.numeric(df1[,1]),as.numeric(df1[,15]), use='pairwise.complete.obs')[[3]]

cor.p<-matrix(NA, nrow=15, ncol=2)
row.names(cor.p)<-colnames(df1)
colnames(cor.p)<-c("cor1","pvalue")

for (i in 1:15)
{
  cor.p[i,1]<-cor.test(as.numeric(df1[,i]),as.numeric(df1[,15]), use='pairwise.complete.obs')[[4]]
  cor.p[i,2]<-cor.test(as.numeric(df1[,i]),as.numeric(df1[,15]), use='pairwise.complete.obs')[[3]]
}

# for ( i in 1:156)

#   
# {plotM<-function(x=1:156, y=157) plot

for (i in 1:15)
{
  png(paste(i,".png", sep=""))
  plot(df1[,i], df1[,15])
  dev.off()
}

#box with jitter
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/miRNA_new_analyses_IBS9232015/wgcna_miRNA")
load("wgcna.part1.files.RData")


colnames(df1)<-gsub("\\|0.007","",colnames(df1))
colnames(df1)<-gsub("\\|0.033","",colnames(df1))
colnames(df1)<-gsub("\\|0","",colnames(df1))
colnames(df1)<-gsub("-","_",colnames(df1))

colnames(df1)

# [1] "hsa_miR_302d_3p" "hsa_miR_631"     "hsa_miR_19a_3p"  "hsa_miR_146b_5p" "hsa_miR_4516"    "hsa_miR_19b_3p" 
# [7] "hsa_miR_378e"    "hsa_miR_1283"    "hsa_miR_324_5p"  "hsa_miR_296_5p"  "hsa_miR_98"      "hsa_miR_1253"   
# [13] "hsa_miR_612"     "hsa_miR_188_5p"  "V2"  

ggplot(na.omit(df1), aes(x=V2,y=hsa_miR_302d_3p))+geom_boxplot(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge", outlier.colour = "red", outlier.shape = 16, outlier.size = 2, notch = F, notchwidth = 0.5)
bar1<-qplot(factor(V2),hsa_miR_302d_3p,data=na.omit(df1), geom=c("boxplot","jitter"), fill=factor(V2), main="BSQ_severity vs hsa-miR-302d-3p",xlab= "BSQ_Severity", ylab="hsa-miR-302d-3p")
x<-bar1 + stat_boxplot(geom ='errorbar')
ggsave(bar1, file="hsa-miR-302d-3p.boxplot.png", width = 8, height = 5, units = "in", dpi = 200)


# remove na from 
df1<-df1[-9,]


for (i in 1:15)
{
  
  x<-qplot(factor(df1[,15]),df1[,i], geom=c("boxplot","jitter"), fill=factor(df1[,15]), main="BSQ_severity",xlab= "BSQ_Severity") 
  y<-x + stat_boxplot(geom ='errorbar')
  print(y)
  ggsave(paste(i,"boxplot.png", sep=""))
}



#########################################################################
microbial dysbiosis index added to mapping file
#import custum index file firmicutes/bacteroidetes ratio

index<-read.delim("custom_index.csv", sep=",", row.names=1)
map3<-merge(map2, md, by="row.names")
row.names(map3)<-map3[,1]
map3<-map3[,c(2:30,32,31)]
t.test(as.numeric(map3$index)~map3$Dx)

md<-read.delim("custom_index.csv", sep=",", row.names=1)
map3<-merge(map2, md, by="row.names")
row.names(map3)<-map3[,1]
map3<-map3[,c(2:30,32,31)]
t.test(as.numeric(map3$index)~map3$Dx)
summary(aov(as.numeric(map3$index)~map3$Bowel.Habit))
summary(aov(as.numeric(map3$index)~map3$cluster))


index<-read.delim("ratios/md.csv", sep=",", row.names=1)
t.test(as.numeric(index$md)~index$Dx)
#data:  as.numeric(index$md) by index$Dx
#t = -0.88145, df = 62.153, p-value = 0.3815
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -15.631074   6.064076
#sample estimates:
# mean in group HC mean in group IBS 
#         38.77419          43.55769 

#test for normality

shapiro.test (ma.phylum[,4])
qqnorm(ma.phylum[,4])
qqline(ma.phylum[,4])
save(exp.clu1.rda, file="R_output/exp.clu1.rda")
save(exp.clu2.rda, file="R_output/exp.clu2.rda")
save(exp.clu3.rda, file="R_output/exp.clu3.rda")
save(exp.hc.rda, file="R_output/exp.hc.rda")

load("exp.clu1.rda")
load("exp.clu2.rda")
load("exp.clu3.rda")
load("exp.hc.rda")
exp.clu1.t<-as.data.frame(t(exp.clu1))
exp.clu2.t<-as.data.frame(t(exp.clu2))
exp.clu3.t<-as.data.frame(t(exp.clu3))
exp.hc.t<-as.data.frame(t(exp.hc))

exp.clu1.t$memb<-c(rep("clu1",9))
exp.clu2.t$memb<-c(rep("clu2",26))
exp.clu3.t$memb<-c(rep("clu3",17))
exp.hc.t$memb<-c(rep("HC",31))

mic.df<-rbind(exp.clu1.t,exp.clu2.t)
wil.12<-t(sapply(mic.df[-72], function(x) 	unlist(wilcox.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
wil.12<-as.data.frame(wil.12)
wil.12$bh<-p.adjust(wil.12$"p.value", method='fdr')
wil.12.sig<-subset(wil.12, p.value<0.05)
save(wil.12.sig, file="wilcox.cluster12.sig.rda")
write.table(wil.12.sig, file="wilcox.cluster12.sig.csv", sep=",", col.names=NA)


mic.df<-rbind(exp.clu2.t,exp.clu3.t)
wil.23<-t(sapply(mic.df[-72], function(x) 	unlist(wilcox.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
wil.23<-as.data.frame(wil.23)
wil.23$bh<-p.adjust(wil.23$"p.value", method='fdr')
wil.23.sig<-subset(wil.23, p.value<0.05)
save(wil.23.sig, file="wilcox.cluster23.sig.rda")
write.table(wil.23.sig, file="wilcox.cluster23.sig.csv", sep=",", col.names=NA)


mic.df<-rbind(exp.clu1.t,exp.clu3.t)
wil.13<-t(sapply(mic.df[-72], function(x) 	unlist(wilcox.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
wil.13<-as.data.frame(wil.13)
wil.13$bh<-p.adjust(wil.13$"p.value", method='fdr')
wil.13.sig<-subset(wil.13, p.value<0.05)
save(wil.13.sig, file="wilcox.cluster13.sig.rda")
write.table(wil.13.sig, file="wilcox.cluster13.sig.csv", sep=",", col.names=NA)


mic.df<-rbind(exp.clu1.t,exp.hc.t)
wil.1.hc<-t(sapply(mic.df[-72], function(x) 	unlist(wilcox.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
wil.1.hc<-as.data.frame(wil.1.hc)
wil.1.hc$bh<-p.adjust(wil.1.hc$"p.value", method='fdr')
wil.1.hc.sig<-subset(wil.1.hc, p.value<0.05)
save(wil.1.hc.sig, file="wilcox.cluster1.hc.sig.rda")
write.table(wil.1.hc.sig, file="wilcox.cluster1.hc.sig.csv", sep=",", col.names=NA)


mic.df<-rbind(exp.clu2.t,exp.hc.t)
wil.2.hc<-t(sapply(mic.df[-72], function(x) 	unlist(wilcox.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
wil.2.hc<-as.data.frame(wil.2.hc)
wil.2.hc$bh<-p.adjust(wil.2.hc$"p.value", method='fdr')
wil.2.hc.sig<-subset(wil.2.hc, p.value<0.05)
save(wil.2.hc.sig, file="wilcox.cluster2.hc.sig.rda")
write.table(wil.2.hc.sig, file="wilcox.cluster2.hc.sig.csv", sep=",", col.names=NA)


mic.df<-rbind(exp.clu3.t,exp.hc.t)
wil.3.hc<-t(sapply(mic.df[-72], function(x) 	unlist(wilcox.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
wil.3.hc<-as.data.frame(wil.3.hc)
wil.3.hc$bh<-p.adjust(wil.3.hc$"p.value", method='fdr')
wil.3.hc.sig<-subset(wil.3.hc, p.value<0.05)
save(wil.3.hc.sig, file="wilcox.cluster3.hc.sig.rda")
write.table(wil.3.hc.sig, file="wilcox.cluster3.hc.sig.csv", sep=",", col.names=NA)



taxa.ibs<-ma.phylum[row.names(ma.phylum)%in%row.names(map2.ibs),]

match(row.names(map2.ibs),row.names(taxa.ibs))
df1<-cbind(taxa.ibs,map2.ibs)
df2<-df1[,-c(72,73,74,75,76,77,101)]
library(ggplot2)
library(reshape2)
df2$Bowel.Habit<-gsub("C", 1, df2$Bowel.Habit)
df2$Bowel.Habit<-gsub("D", 2, df2$Bowel.Habit)
df2$Bowel.Habit<-gsub("M", 6, df2$Bowel.Habit)
df2$Bowel.Habit<-as.numeric(df2$Bowel.Habit)
df2$Gender<-gsub("M", 1, df2$Gender)
df2$Gender<-gsub("F", 2, df2$Gender)
df2$Gender<-as.numeric(df2$Gender)
df2$Dx<-NULL
for ( i in 72:93){
df2[,i]<-as.numeric(as.character(df2[,i]))
}

ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
lm.D90 <- lm(weight ~ group - 1) # omitting intercept

anova(lm.D9)
summary(lm.D90)

mic.df<-df2[,c(1:71,72)]

lm.bh<-summary(lm(df2[,1] ~ df2[,72]))$coefficients[,4][2]
lm.bh<-t(sapply(df2[,-72], function(x) summary(lm(x~df2$Bowel.Habit)$coefficients[,4][2])))
p.lm<-matrix(NA, nrow=71, ncol=21)

for ( i in 1:71)

{
p.lm[i,1]<-summary(lm(df2[,i]~df2[,72]))$coefficients[,4][2]
p.lm[i,2]<-summary(lm(df2[,i]~df2[,73]))$coefficients[,4][2]
p.lm[i,3]<-summary(lm(df2[,i]~df2[,74]))$coefficients[,4][2]
p.lm[i,4]<-summary(lm(df2[,i]~df2[,75]))$coefficients[,4][2]
p.lm[i,5]<-summary(lm(df2[,i]~df2[,76]))$coefficients[,4][2]
p.lm[i,6]<-summary(lm(df2[,i]~df2[,77]))$coefficients[,4][2]
p.lm[i,7]<-summary(lm(df2[,i]~df2[,78]))$coefficients[,4][2]
p.lm[i,8]<-summary(lm(df2[,i]~df2[,79]))$coefficients[,4][2]
p.lm[i,9]<-summary(lm(df2[,i]~df2[,80]))$coefficients[,4][2]
p.lm[i,10]<-summary(lm(df2[,i]~df2[,81]))$coefficients[,4][2]
p.lm[i,11]<-summary(lm(df2[,i]~df2[,82]))$coefficients[,4][2]
p.lm[i,12]<-summary(lm(df2[,i]~df2[,83]))$coefficients[,4][2]
p.lm[i,13]<-summary(lm(df2[,i]~df2[,84]))$coefficients[,4][2]
p.lm[i,14]<-summary(lm(df2[,i]~df2[,85]))$coefficients[,4][2]
p.lm[i,15]<-summary(lm(df2[,i]~df2[,86]))$coefficients[,4][2]
p.lm[i,16]<-summary(lm(df2[,i]~df2[,87]))$coefficients[,4][2]
p.lm[i,17]<-summary(lm(df2[,i]~df2[,88]))$coefficients[,4][2]
p.lm[i,18]<-summary(lm(df2[,i]~df2[,89]))$coefficients[,4][2]
p.lm[i,19]<-summary(lm(df2[,i]~df2[,90]))$coefficients[,4][2]
p.lm[i,20]<-summary(lm(df2[,i]~df2[,91]))$coefficients[,4][2]
p.lm[i,21]<-summary(lm(df2[,i]~df2[,92]))$coefficients[,4][2]
}

colnames(p.lm)<-colnames(df2[,72:92])
row.names(p.lm)<-colnames(df2[,1:71])

write.table(p.lm, file="R_output/p.lm.traits.csv")
o.sig<-subset(p.lm, p.lm[,18]<0.05)

wil.3.hc<-as.data.frame(wil.3.hc)
wil.3.hc$bh<-p.adjust(wil.3.hc$"p.value", method='fdr')
wil.3.hc.sig<-subset(wil.3.hc, p.value<0.05)
save(wil.3.hc.sig, file="wilcox.cluster3.hc.sig.rda")
write.table(wil.3.hc.sig, file="wilcox.cluster3.hc.sig.csv", sep=",", col.names=NA)


qplot(x=Var1, y=Var2, data=melt(cor(df2,use = 'pairwise.complete.obs')), fill=value, geom="tile")
                                   
cor.test(df2$Bowel.Habit,df2$'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium; s__.1', use = 'pairwise.complete.obs')
                                    
                                 
boxplot(df2$Bowel.Habit,df2$'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__; s__')
                                   
phylum.ibs.pre<-as.data.frame(phylum.ibs.pre)
                                   
phylum.ibs.pre.fbratio<-phylum.ibs.pre$'k__Bacteria; p__Firmicutes'/phylum.ibs.pre$'k__Bacteria; p__Bacteroidetes' 
mean(phylum.ibs.pre.fbratio)
                               # 0.7482137
                               
                               ibs.norLike1<-as.data.frame(t(ibs.norLike))
                               
phylum.ibs.norlike.fbratio<-ibs.norLike1$'k__Bacteria; p__Firmicutes'/ibs.norLike1$'k__Bacteria; p__Bacteroidetes' 
                               
mean(phylum.ibs.norlike.fbratio)
                               # 0.8671705
                               
exp.hc<-as.data.frame(exp.hc)
phylum.exp.hc.fbratio<-exp.hc$'k__Bacteria; p__Firmicutes'/exp.hc$'k__Bacteria; p__Bacteroidetes' 
                               
mean(phylum.exp.hc.fbratio)
# 0.9491604
                                   
                               
                # new mapping file
                               hc1<-subset(map2, map2$Dx=="HC")
                               clust1=c(rep("hc",31)
                               memb2<-rbind(memb1,cbind(row.names(hc1),clust1))
                               memb2$clust1<-gsub("1","IBS_like",memb2$clust1)
                               memb2$clust1<-gsub("2","Normal_like",memb2$clust1)
                               row.names(memb2)<-memb2[,1]
                               memb2<-memb2[row.names(map2),]
                               match(row.names(map2), memb2[,1])
                               map3<-cbind(map2,memb2[,2])
                               colnames(map3)[44]<-"Clust1"
                               map3<-map3[,c(1:42,44,43)]
                          write.table(map3, file="swpna.cluster.mapping.txt", sep='\t')     
                          
#############################
map4<- subset( map3,  map3$Clust1=="IBS_like"| map3$Clust1=="Normal_like")
map5<-map4[,-c(1,2,3,4,5,6,7,10,12:20,34:42,44)]
                          
ttest.dis.cont<-t(sapply(map5[,-17], function(x) 	unlist(t.test(x~map5$Clust1, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
# ttest.dis.cont$logFC<-ttest.dis.cont$'estimate.mean in group Relapse'-ttest.dis.cont$'estimate.mean in group No Relapse'
# ttest.dis.cont$abslogFC<-abs(ttest.dis.cont$logFC)
write.table( ttest.dis.cont, file="t-test-pheno-groupsibs.csv", col.names=NA, sep=",")
                         
                         # ttest.dis.cont$absFC<-2^ttest.dis.cont$abslogFC
                         ttest.dis.cont$FC<-NA
                         
                         for (i in 1: 7) {
                           if(ttest.dis.cont[i,2]/ttest.dis.cont[i,1]>1)
                           {ttest.dis.cont[i,8]<-ttest.dis.cont[i,2]/ttest.dis.cont[i,1]}
                           else if (ttest.dis.cont[i,2]/ttest.dis.cont[i,1]<=1)
                           {ttest.dis.cont[i,8]<--1*(ttest.dis.cont[i,1]/ttest.dis.cont[i,2])}
                         }
                         
                         save(ttest.dis.cont, file="ttest.norlikeibs.vs.controls.rda")
                         write.table(ttest.dis.cont, file="ttest.norlikeibs.vs.controls.csv", sep=",", col.names=NA)
                         ibs.pre.sig<-ibs.pre[,colnames(ibs.pre)%in%row.names(ttest.dis.cont.sig)]
                         ibs.norLike.sig<-ibs.norLike[,colnames(ibs.norLike)%in%row.names(ttest.dis.cont.sig)]
                         exp.hc.sig<-exp.hc[,colnames(exp.hc)%in%row.names(ttest.dis.cont.sig)]
                         
    ##################################################################################
##phyloseq:
                         install.packages("biom")
                         library(biom)
header1<-scan("otu_table_SCOR_biopsy.txt", nlines = 2, what = character())[8:91]
genus=read.table("otu_table_SCOR_biopsy.txt", sep="\t", skip=2, header=FALSE, row.names = 1)	
names(genus)<-header1    
taxa.names = genus$ConsensusLineage

taxa.mat<-matrix(c(unlist(strsplit(as.character(genus$ConsensusLineage),";"))), ncol=7, byrow=TRUE)
row.names(taxa.mat)<-row.names(genus)
colnames(taxa.mat)<-c("kingdom","Phylum","class","order","family","genus","species")
save(taxa.mat, file="taxa.mat.rda")
load("ma.phylum.rda")
microbMat<-ma.phylum

library(BiocInstaller)
biocLite("phyloseq")
filepath = system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", package="phyloseq")
library("phyloseq"); packageVersion("phyloseq")
ibs.biom<-microbio_me_qiime(genus)
import_qiime(genus)
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/IBSmicrobiomeProject/biopsy_scor_microbiome/jonathan_analysis/")
otufile=("otu_table_SCOR_biopsy.biom") 
file <-read_biom("otu_table_SCOR_biopsy.biom")      

##### integration

# import miRNA data:
  # import microbiome data
# make the matrices with same samples and same order

load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/TissueRCC/norm.data4.rda")
mirna<-norm.data4
colnames(mirna)<-substr(colnames(mirna),1,5)
microbMat<-df2[,1:71]
row.names(microbMat)<-substr(row.names(microbMat),8,12)
microbMat1<-t(microbMat)
microbMat2<-microbMat1[,colnames(microbMat1)%in%colnames(mirna)]
mirna1<-mirna[,colnames(mirna)%in%colnames(microbMat2)]
microbMat2<-microbMat2[,colnames(mirna1)]
match(colnames(microbMat2),colnames(mirna1)) #only IBS will have 1:29
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
# [38] 38 39 40 41 42 43
library(mixOmics)
color<-color.jet
par(mfrow = c(1,2));
png("R_output/mixOmics.png", height=2000, width=3000, res=300)
liver.splsda <- imgCor(t(mirna1), t(microbMat2), X.var.names = TRUE,
                       Y.var.names = TRUE,
                       sideColors = TRUE,
                       interactive.dev = TRUE,
                       main = TRUE)
                       # color, row.cex, col.cex,symkey, keysize,
                       # xlab, ylab, margins, lhei, lwid)
dev.off()

data1<-cbind(t(mirna1),t(microbMat2))
library(ellipse)
ctab <- cor(data1)
corMat<-round(ctab, 2)
# sub1<-function(x)(subset(x, x[i]>=0.5))
highCor<-matrix(NA, nrow=561, ncol=561)
row.names(highCor)<-row.names(corMat)
colnames(highCor)<-colnames(corMat)
for ( i in 1:dim(corMat)[1])
  for ( j in 1:dim(corMat)[2])
  {
    if(corMat[i,j]>0.6)
    {highCor[i,j]<- 1}
    else if(corMat[i,j]< -0.6)
    {highCor[i,j]<- 1}
    else (highCor[i,j]<- 0)
  }

mirnaMicroCor<-highCor[1:490,491:561]
mirnaMicroCor<-as.data.frame(mirnaMicroCor)

mirnaMicroCor$sum1<-apply(mirnaMicroCor,1,sum)

mirnaMicroCor1<-mirnaMicroCor[mirnaMicroCor$sum1>0,]
data1<-as.data.frame(data1)


plot(data1$'hsa-miR-196a-5p|0', data1$'k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__  AND k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__.1')

plot(data1$'hsa-miR-574-5p|0', data1$'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__')

plot(data1$'hsa-miR-130b-3p|0', data1$'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__')



write.table(mirnaMicroCor1, file="R_output/mirnaMicroCor1.csv", sep=",", col.names=NA)
plot(data1$'hsa-miR-34c-5p|0', data1$'k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Sutterella')


#difference in means of clusters

aov(map2$cluster, map2$PSS)
