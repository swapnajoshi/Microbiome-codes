####################
# author: Swapna 
# date  : 4/25/2016
####################

setwd("C:/Users/swapnajoshi-admin/Documents/swapna") 

#WGCNA
library(WGCNA)
library(flashClust)
library(vegan)

#Heatmap
library(dplyr)
library(NMF)
library(RColorBrewer)

# Import distance metric
unwei.dm <- read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/IBSmicrobiomeProject/IBSmicrobiomeAnalysisR/data/rawData/dm/unweighted_unifrac_dm.txt", row.names=1)
# import mapping file
map1 <- read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/IBSmicrobiomeProject/IBSmicrobiomeAnalysisR/data/rawData/ibs.mapping_new.csv",  sep=',', row.names=1)
map2<-map1[row.names(map1)%in%colnames(unwei.dm),] ## 2 samples did not have data
map2<-map2[colnames(unwei.dm),]
match(row.names(map2),colnames(unwei.dm))
unwei.dm1<-as.dist(unwei.dm)
sampleTree = flashClust(unwei.dm1, method = "complete")
labels1=map2$Dx
dat.ordered<-unwei.dm[sampleTree$order,]

sizeGrWindow(15,12)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Heirarchical Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
abline(h = 0.6, col = "red")

clust1 = cutreeStatic(sampleTree, cutHeight = 0.75, minSize = 1)
table(clust1)
lab<-map2[sampleTree$labels,]$Dx
lab1<-map2[sampleTree$labels,]$NDPNum
df2<-as.data.frame(cbind(clust1,as.character(lab)))
head(df2)
set.seed(123)
rainbow(n=6)
# [1] "#FF0000FF" "#FFFF00FF" "#00FF00FF" "#00FFFFFF" "#0000FFFF" "#FF00FFFF"

df2$color1<-NA
for ( i in 1: 83)
{if (df2[i,1]=="1")
{df2[i,3]<-"red"}
  else if (df2[i,1]=="2")
  {df2[i,3]<-"cyan"}
  else if (df2[i,1]=="3")
  {df2[i,3]<-"orange"}
}

df2$colorDx<-NA
for ( i in 1: 83)
{if (df2[i,2]=="IBS")
{df2[i,4]<-"cyan"}
  else if (df2[i,2]=="HC")
  {df2[i,4]<-"red"}
}
df2<-cbind(df2, as.character(lab1))

library(dendextend)
png("unweightedUnifracDistCompleteHierarchical.png", res = 200, height = 2000, width = 2500)
plot(sampleTree, main = "Heirachical Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
abline(h = 0.75, col = "red")
colored_bars(colors= df2[,c(3,4)], dend = sampleTree, y_shift = 0.4, rowLabels = c("Cluster", "Dx"),sort_by_labels_order 	
= TRUE)
dev.off()

# Cluster membership with IBS as well as healthy controls

mapIbsHcCluster<-map2
mapIbsHcCluster$IbsHcClusterMemb<-clust1

write.table(mapIbsHcCluster, file="mapIbsHcClusterBig.csv",sep=",", col.names=NA)
colnames(mapIbsHcCluster)
littleIbsMappingIBsHcClust<-mapIbsHcCluster[,c(10,11,15,16,17,18,20,30:51,76:79)]
write.table(littleIbsMappingIBsHcClust, file="littleIbsMappingIBsHcClust.csv",sep=",", col.names=NA)

save(df2, file="df2.rda")

####################################################

# Are there clusters within IBS? 

# Import distance metric from IBS only bidiv

unwei.dm <- read.delim("C:/Users/swapnajoshi-admin/Documents/swapna/qiime_output/IBS_bdiv_even_9188/unweighted_unifrac_dm.txt", row.names=1)
# import mapping file
map1 <- read.delim("C:/Users/swapnajoshi-admin/Documents/swapna/ibs.mapping_new.csv",  sep=',', row.names=1)
map2<-map1[row.names(map1)%in%colnames(unwei.dm),] ## 2 samples did not have data
map2<-map2[colnames(unwei.dm),]
match(row.names(map2),colnames(unwei.dm))
unwei.dm1<-as.dist(unwei.dm)
sampleTree = flashClust(unwei.dm1, method = "ward")
labels1=map2$Dx
dat.ordered<-unwei.dm[sampleTree$order,]

sizeGrWindow(15,12)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Heirarchical Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
abline(h = 1.2, col = "red")

clust1 = cutreeStatic(sampleTree, cutHeight = 1.2, minSize = 1)
table(clust1)
lab<-map2[sampleTree$labels,]$Dx
lab1<-map2[sampleTree$labels,]$NDPNum
df2<-as.data.frame(cbind(clust1,as.character(lab)))
head(df2)
set.seed(123)
rainbow(n=6)
# [1] "#FF0000FF" "#FFFF00FF" "#00FF00FF" "#00FFFFFF" "#0000FFFF" "#FF00FFFF"

df2$color1<-NA
for ( i in 1: 52)
{if (df2[i,1]=="1")
{df2[i,3]<-"#FF0000FF"}
  else if (df2[i,1]=="2")
  {df2[i,3]<-"#00FF00FF"}
  else if (df2[i,1]=="3")
  {df2[i,3]<-"#0000FFFF"}
}

df2$colorDx<-NA
for ( i in 1: 52)
{if (df2[i,2]=="IBS")
{df2[i,4]<-"#0066FFFF"}
  else if (df2[i,2]=="HC")
  {df2[i,4]<-"#00FFFFFF"}
}
df2<-cbind(df2, as.character(lab1))

library(dendextend)
png("unweightedUnifracDistCompleteHierarchicalIBS.png", res = 300, height = 2000, width = 2500)
plot(sampleTree, main = "Heirachical Clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
abline(h = 1.2, col = "red")
colored_bars(colors= df2[,c(3,4)], dend = sampleTree, y_shift = 0.3, rowLabels = c("Cluster", "Dx"),sort_by_labels_order 	
             = TRUE)
dev.off()

# Cluster membership within IBS

mapIbsCluster<-map2
mapIbsCluster$IbsClusterMemb<-clust1

write.table(mapIbsCluster, file="mapIbsClusterBig.csv",sep=",", col.names=NA)
colnames(mapIbsCluster)
littleIbsMappingIBsClust<-mapIbsCluster[,c(10,11,15,16,17,18,20,30:51,76:79)]
write.table(littleIbsMappingIBsClust, file="littleIbsMappingIBsClust.csv",sep=",", col.names=NA)

save(df2, file="df2.Ibs.rda")

# found 2 clusters

################################################################
# multi-omics
################################################################

library(vegan)

setwd("C:/Users/swapnajoshi-admin/Documents/swapna") 
# read OTU table
header1<-scan("otu_table_even_9188.biom.txt", nlines = 2, what = character())[8:91]
genus=read.table("otu_table_even_9188.biom.txt", sep="\t", skip=2, header=FALSE, row.names = 1)	
names(genus)<-header1
genus1 = as.matrix(genus[,-dim(genus)[2]])
genus2 = scale(genus1, center=F, scale=colSums(genus1)) ## sacle it so they are relative abundances and not counts
genus2 = as.data.frame(t(genus2))
save(genus2, file="genus2.rda")
taxa.names = genus$ConsensusLineage

save(genus1, file="genus1.rda")
save(genus, file="genus.rda")

extract.name.level = function(x, level){
  a=c(unlist(strsplit(x,';')),'Other')
  paste(a[1:min(level,length(a))],collapse=';')
}

otu2taxonomy = function(x, level, taxa=NULL){
  if(is.null(taxa))
  {taxa = colnames(x)}
  if(length(taxa)!=dim(x)[2])
  {print("ERROR: taxonomy should have the same length
         as the number of columns in OTU table")
    return;}
  level.names = sapply(as.character(taxa),
                       function(x)
                         extract.name.level(x,level=level))
  t(apply(x, 1, function(y) tapply(y,level.names,sum)))
}

d.species = otu2taxonomy(genus2,level=7,taxa=taxa.names)
dim(d.species)
d.genus = otu2taxonomy(genus2,level=6,taxa=taxa.names)
dim(d.genus)
d.family = otu2taxonomy(genus2,level=5,taxa=taxa.names)
dim(d.family)
d.order = otu2taxonomy(genus2,level=4,taxa=taxa.names)
dim(d.order)
d.class = otu2taxonomy(genus2,level=3,taxa=taxa.names)
dim(d.class)
d.phylum = otu2taxonomy(genus2,level=2,taxa=taxa.names)
dim(d.phylum)

#########################################################

dim(d.species)

# filtering of probes can be either based on means or standard deviations
# based on means
ma.species = subset(t(d.species), rowMeans(t(d.species)) > 0.01); dim(ma.species)

# based on standard deviation
dt.species<-as.data.frame(t(d.species))
dt.species$sd1<-apply(dt.species,1, sd)
hist(dt.species$sd1)

sd.species<-subset(dt.species, dt.species$sd1>=0.02); dim(sd.species)
sd.species<-as.data.frame(sd.species)
sd.species$sd1<-NULL

datExpr0<-t(ma.species)
save(ma.species, file="ma.species.rda")

load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/TissueRCC/norm.data4.rda"); dim(norm.data4)
mirna<-norm.data4
colnames(mirna)<-substr(colnames(mirna),1,5)
microbMat<-datExpr0
row.names(microbMat)<-substr(row.names(microbMat),8,12)
microbMat1<-t(microbMat)
microbMat2<-microbMat1[,colnames(microbMat1)%in%colnames(mirna)]
mirna1<-mirna[,colnames(mirna)%in%colnames(microbMat2)]
microbMat2<-microbMat2[,colnames(mirna1)]
match(colnames(microbMat2),colnames(mirna1))
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
# [38] 38 39 40 41 42 43 44
library(mixOmics)
color<-color.jet
par(mfrow = c(1,2));
png("mixOmics.png", height=2000, width=3000, res=300)
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
ctab <- cor(data1, method = "spearman")

corMat<-round(ctab, 2)
# sub1<-function(x)(subset(x, x[i]>=0.5))
highCor<-matrix(NA, nrow=356, ncol=356)
row.names(highCor)<-row.names(corMat)
colnames(highCor)<-colnames(corMat)
for ( i in 1:dim(corMat)[1])
  for ( j in 1:dim(corMat)[2])
  {
    if(corMat[i,j] >= 0.4){
    highCor[i,j] <- 1
    }   else if(corMat[i,j] < -0.4) {
    highCor[i,j] <- 1
    } else {
      highCor[i,j] <- 0
    }
  }

mirnaMicroCor<-highCor[1:339,340:356]
mirnaMicroCor<-as.data.frame(mirnaMicroCor)

mirnaMicroCor$sum1<-apply(mirnaMicroCor,1,sum)
mirnaMicroCor1<-mirnaMicroCor[mirnaMicroCor$sum1>0,]
selSp<-apply(mirnaMicroCor1,2,sum)
selSp1<-names(selSp[selSp>0])
mirnaMicroCor2 <- mirnaMicroCor1[,selSp1]

dim(mirnaMicroCor2)

write.table(mirnaMicroCor2, file="mirnaMicroCor2.csv", sep=",", col.names=NA)

bacSel <- colnames(mirnaMicroCor1)
mirSel <- row.names(mirnaMicroCor1)

cor.microb <- t(microbMat2[row.names(microbMat2)%in%bacSel,]); dim(cor.microb)
cor.mirna<-t(mirna1[row.names(mirna1) %in% mirSel,]); dim(cor.mirna)

data2<-as.data.frame(cbind(cor.microb,cor.mirna))
littleIbsMappingIBsHcClust<-read.delim(file="littleIbsMappingIBsHcClust.csv",sep=",", row.names=1); dim(littleIbsMappingIBsHcClust)
map.mirna<-littleIbsMappingIBsHcClust[substr(row.names(littleIbsMappingIBsHcClust),8,13)%in%row.names(data2),]; dim(map.mirna)

row.names(map.mirna)<-substr(row.names(map.mirna),8,13)
map.mirna<-map.mirna[row.names(data2),]
match(row.names(data2),row.names(map.mirna))
data3<-cbind(data2, map.mirna[,c(5,3,4,33)])
colnames(data3)<-gsub(" ","",colnames(data3))
colnames(data3)<-gsub(";","",colnames(data3))
colnames(data3)<-gsub("-","_",colnames(data3))
colnames(data3)<-gsub("\\|","",colnames(data3))
colnames1<-colnames(data3)

setwd("C:/Users/swapnajoshi-admin/Documents/swapna/microbMiRNAcorPlots/")
for ( i in 1:17)
  for (j in 18:68) {
    cat(paste("plot",i,j))
    P<-ggplot(data3, aes(data3[,i],data3[,j])) + geom_point(size=4) 
    q<-P+geom_point(aes(color=as.factor(data3$Group)))+ xlab(colnames1[i]) + ylab(colnames1[j])
ggsave(q, filename = paste(colnames1[i],colnames1[j],".png"))
  }

P<-ggplot(data3, aes(data3[,9],data3[,38])) + geom_point(size=4) 
q<-P+geom_point(aes(color=as.factor(data3$Group)))+ xlab("F_Lachnospiraceae") + ylab("hsa_miR_200c_3p0") + geom_smooth(aes(data3[,9],data3[,38],colour=factor(Group)), method=lm, se=FALSE)
q
ggsave(q, filename = "hsa_miR_200c_3p0_F_Lachnospiraceae_lm.png")

P<-ggplot(data3, aes(data3[,7],data3[,19])) + geom_point(size=4) 
q<-P+geom_point(aes(color=as.factor(data3$Group)))+ xlab("F_Alicyclobacillus") + ylab("hsa_miR_195_5p0") + geom_smooth(aes(data3[,7],data3[,19],colour=factor(Group)), method=lm, se=FALSE)
q
ggsave(q, filename = "hsa_miR_195_5p0_F_Alicyclobacillus_lm.png")

################################

# pubmed search for articles related to the miRNAs od interest
# hsa-miR-141-3p|0
# hsa-miR-144-3p|0.005
# hsa-miR-125b-5p|0.009
# hsa-miR-139-5p|0
# hsa-miR-204-5p|0
# hsa-miR-100-5p|0
# hsa-miR-1|0
# hsa-let-7g-5p|0
# hsa-miR-411-5p|0
# hsa-miR-200c-3p|0
# hsa-miR-374a-5p|0

library(RISmed)
res <- EUtilsSummary('miR-200c, microbiome', type='esearch', db='pubmed')
summary(res)

# Query:
#   miR-200c[All Fields] AND ("microbiota"[MeSH Terms] OR "microbiota"[All Fields] OR "microbiome"[All Fields]) 
# 
# Result count:  1

###############################

setwd("C:/Users/swapnajoshi-admin/Documents/swapna") 

littleIbsMappingIBsClust$Description<-NULL
littleIbsMappingIBsClust$Group<-NULL
littleIbsMappingIBsClust$Gender<-gsub("M","1",littleIbsMappingIBsClust$Gender)
littleIbsMappingIBsClust$Gender<-gsub("F","2",littleIbsMappingIBsClust$Gender)
map4<-apply(littleIbsMappingIBsClust,2,as.numeric)
row.names(map4)<-row.names(littleIbsMappingIBsClust)

save(map4, file="map.numeric.clu1.2.rda")

p.traits1<-matrix(NA, nrow=31, ncol=1)
for ( i in 1: 30)
{p.traits1[i,1]<-summary(lm(map4[,i] ~ map4[,31]))$coefficients[,4][2]}
row.names(p.traits1)<-colnames(map4)
p.traits1
write.table(p.traits1, file="TraitsIbsClusters_lm.csv", col.names=NA, sep=",")

boxplot(map4[,9] ~ as.factor(map4[,31]))

# ##################### firmicutes bacteroidetes ratio at phylum level
# 
# ma.phylum1<-t(ma.phylum)
# match(row.names(ma.phylum1), row.names(littleIbsMapping))
# littleIbsMapping<-littleIbsMapping[row.names(ma.phylum1),]
# ma.phylum1<-as.data.frame(ma.phylum1)
# littleIbsMapping$FBRatio<-ma.phylum1$'k__Bacteria; p__Firmicutes'/ma.phylum1$'k__Bacteria; p__Bacteroidetes'
# littleIbsMapping$FBRatio<-as.numeric(as.character(littleIbsMapping$FBRatio))
# littleIbsMapping<-as.data.frame(littleIbsMapping)
# 
# t.test(littleIbsMapping[,31]~ littleIbsMapping[,7])
# boxplot(littleIbsMapping[,31]~ littleIbsMapping[,7])
# 
# #mean in group HC mean in group IBS 
# #       0.9491604         0.7893911  p-value = 0.4082
# 
# 
# summary(aov(littleIbsMapping[,31]~ littleIbsMapping[,8]))
# littleIbsMapping[,8]<-gsub("U","M",littleIbsMapping[,8])
# #p=0.807
# boxplot(littleIbsMapping[,31]~ littleIbsMapping[,8])
# summary(aov(littleIbsMapping[,31]~ littleIbsMapping[,25]))
# #7.17e-07
# boxplot(littleIbsMapping[,31]~ littleIbsMapping[,25])
# 
# 
# boxplot(map4[,31]~ map4[,7])
# 
# t.test(littleIbsMapping[,31]~ littleIbsMapping[,9])
# 
# 
# 
# t.test(littleIbsMapping[littleIbsMapping$Dx=="IBS",31]~ littleIbsMapping[littleIbsMapping$Dx=="IBS",9])
# boxplot(littleIbsMapping[littleIbsMapping$Dx=="IBS",31]~ littleIbsMapping[littleIbsMapping$Dx=="IBS",9])
# 
# t.test(littleIbsMapping[littleIbsMapping$Dx=="HC",31]~ littleIbsMapping[littleIbsMapping$Dx=="HC",9])
# boxplot(littleIbsMapping[littleIbsMapping$Dx=="HC",31]~ littleIbsMapping[littleIbsMapping$Dx=="HC",9])
# 
# 

#######################################
# alicyclobacillus membership

setwd("C:/Users/swapnajoshi-admin/Documents/swapna") 

match(row.names(datExpr0), row.names(littleIbsMappingIBsHcClust))
littleIbsMappingIBsHcClust$Description<-NULL
# littleIbsMappingIBsHcClust$Group<-NULL
littleIbsMappingIBsHcClust$Gender<-gsub("M","1",littleIbsMappingIBsClust$Gender)
littleIbsMappingIBsHcClust$Gender<-gsub("F","2",littleIbsMappingIBsClust$Gender)
map4<-apply(littleIbsMappingIBsHcClust,2,as.numeric)
row.names(map4)<-row.names(littleIbsMappingIBsHcClust)
map4<-as.data.frame(map4)
map4$Alicyclobacillus<-NA
for (i in 1:83) {
  if(datExpr0[i,7] >= 0.005) {
    map4[i,33] <- 1
  } else if(datExpr0[i,7] < 0.005) {
    map4[i,33] <- 0
  } 
}
table(map4$Alicyclobacillus)

map4.ibs<-subset(map4, map4$Group == 2)
save(map4.ibs, file = "IBSmappingAlicyclobacillus.rda")
save(map4, file = "mappingAlicyclobacillus.rda")

map4.ibs$Group<-NULL
ab.traits1<-matrix(NA, nrow=32, ncol=1)
for ( i in 1: 30)
{ab.traits1[i,1]<-summary(lm(map4.ibs[,i] ~ map4.ibs[,32]))$coefficients[,4][2]}
row.names(ab.traits1)<-colnames(map4.ibs)
ab.traits1
write.table(ab.traits1, file="alicyclobacillusTraitsIbs_lm.csv", col.names=NA, sep=",")

boxplot(map4.ibs[,9] ~ as.factor(map4.ibs[,32]))


##############################################trait correlation with membership ibs and hc

map4<-read.delim("C:/Users/swapnajoshi-admin/Documents/swapna/littleIbsMappingIBsHcClust.csv", sep=",", row.names = 1)
for(i in 1:83) { if(map4[i,3] == 1) {map4[i,33] <- 4}}

map4$Description<-NULL
map4$Gender<-gsub("M","1",map4$Gender)
map4$Gender<-gsub("F","2",map4$Gender)

# between cluster 1 and 2

map4_1<-subset(map4, map4[,32] == 1 | map4[,32] == 2)
traits1<-matrix(NA, nrow=32, ncol=1)
for ( i in 1: 31)
{traits1[i,1]<-summary(lm(map4_1[,i] ~ map4_1[,32]))$coefficients[,4][2]}
row.names(traits1)<-colnames(map4_1)
traits1
# 
# [,1]
# Gender                            0.65818938
# Age_orig                          0.41016023
# Group                             0.38648787
# BH                                0.28960571
# Sex                               0.65818938
# Age                               0.41798435
# BMI                               0.59786452
# BSQ_OverallSx                     0.49765472
# BSQ_AbdPain                       0.52240501
# BSQ_Bloating                      0.03165633
# BSQ_UsualSeverity                 0.01949026
# BSQ_AgeOnset                      0.37077189
# BSQdrv_SxDurationYrs              0.58500187
# ETI_General_Score                 0.41930571
# ETI_Physical_Score                0.27578707
# ETI_Emotional_Score               0.44812027
# ETI_Sexual_Score                  0.40150547
# ETI_Total_Score                   0.29311518
# HAD_Anxiety                       0.45433811
# HAD_Depression                    0.10957752
# ACE_Emotional_Abuse               0.86212364
# ACE_Physical_Abuse                0.73673763
# ACE_Sexual_Abuse                  0.76044454
# ACE_Substance_Abuse               0.55040365
# ACE_Parental_DivorceSep           0.83926435
# ACE_Household_Mental_Illness      0.69772954
# ACE_Incarcerated_Household_Member 0.38828990
# ACE_Parents_Treated_Violently     0.37085757
# ACE_Score                         0.67131358
# VSI_Score                         0.00783894
# PSS_Score                         0.72753955
# IbsHcClusterMemb                          NA
write.table(traits1, file="ibsHcCluster12Traits_lm.csv", col.names=NA, sep=",")

p <- ggplot(map4, aes(factor(IbsHcClusterMemb), BSQ_Bloating))
q <- p + geom_boxplot(aes(fill = factor(IbsHcClusterMemb)))
ggsave(q, file="BSQ_Bloating_ibsHcCluster.png")

p <- ggplot(map4, aes(factor(IbsHcClusterMemb), BSQ_UsualSeverity))
q <- p + geom_boxplot(aes(fill = factor(IbsHcClusterMemb)))
ggsave(q, file="BSQ_UsualSeverity_ibsHcCluster.png") ## missing points

p <- ggplot(map4, aes(factor(IbsHcClusterMemb), VSI_Score))
q <- p + geom_boxplot(aes(fill = factor(IbsHcClusterMemb)))
ggsave(q, file="VSI_Score_ibsHcCluster.png")


p <- ggplot(map4, aes(factor(IbsHcClusterMemb), BSQ_OverallSx))
q <- p + geom_boxplot(aes(fill = factor(IbsHcClusterMemb)))
ggsave(q, file="BSQ_OverallSx_ibsHcCluster.png") ### missing points

# between cluster 2 and 3 

map4_2<-subset(map4, map4[,32] == 2 | map4[,32] == 3);dim(map4_2)
traits2<-matrix(NA, nrow=32, ncol=1)
for ( i in 1: 31)
{traits2[i,1]<-summary(lm(map4_2[,i] ~ map4_2[,32]))$coefficients[,4][2]}
row.names(traits2)<-colnames(map4_2)
traits2
# Gender                            0.61704237
# Age_orig                          0.84625627
# Group                             0.39125583
# BH                                0.27351403
# Sex                               0.61704237
# Age                               0.85678589
# BMI                               0.35616743
# BSQ_OverallSx                     0.54711595
# BSQ_AbdPain                       0.56294021
# BSQ_Bloating                      0.04092474
# BSQ_UsualSeverity                 0.44214311
# BSQ_AgeOnset                      0.22098606
# BSQdrv_SxDurationYrs              0.49429431
# ETI_General_Score                 0.82715094
# ETI_Physical_Score                0.18552862
# ETI_Emotional_Score               0.94625095
# ETI_Sexual_Score                  0.24549289
# ETI_Total_Score                   0.48013156
# HAD_Anxiety                       0.60475465
# HAD_Depression                    0.26017662
# ACE_Emotional_Abuse               0.67618761
# ACE_Physical_Abuse                0.23677723
# ACE_Sexual_Abuse                  0.41870262
# ACE_Substance_Abuse               0.61587644
# ACE_Parental_DivorceSep           0.15162514
# ACE_Household_Mental_Illness      0.68239615
# ACE_Incarcerated_Household_Member 0.49797988
# ACE_Parents_Treated_Violently     0.30177469
# ACE_Score                         0.33245359
# VSI_Score                         0.23154178
# PSS_Score                         0.82325559
# IbsHcClusterMemb                          NA

write.table(traits2, file="ibsHcCluster23Traits_lm.csv", col.names=NA, sep=",")

# between cluster 1 and 3 

map4_3<-subset(map4, map4[,32] == 1 | map4[,32] == 3);dim(map4_3)
traits3<-matrix(NA, nrow=32, ncol=1)
for ( i in 1: 31)
{traits3[i,1]<-summary(lm(map4_3[,i] ~ map4_3[,32] + map4_3[,6]))$coefficients[,4][2]}
row.names(traits3)<-colnames(map4_3)
traits3
# Gender                            0.67864689
# Age_orig                          0.14025897
# Group                             0.43282517
# BH                                0.82261400
# Sex                               0.67864689
# Age                               1.00000000
# BMI                               0.44197888
# BSQ_OverallSx                     0.04780228
# BSQ_AbdPain                       0.12285288
# BSQ_Bloating                      0.76155543
# BSQ_UsualSeverity                 0.07255482
# BSQ_AgeOnset                      0.10202509
# BSQdrv_SxDurationYrs              0.59436306
# ETI_General_Score                 0.54624582
# ETI_Physical_Score                0.92168128
# ETI_Emotional_Score               0.32887446
# ETI_Sexual_Score                  0.50766027
# ETI_Total_Score                   0.61866226
# HAD_Anxiety                       0.23284959
# HAD_Depression                    0.69075972
# ACE_Emotional_Abuse               0.83063897
# ACE_Physical_Abuse                0.51376430
# ACE_Sexual_Abuse                  0.40037562
# ACE_Substance_Abuse               0.77112047
# ACE_Parental_DivorceSep           0.25892730
# ACE_Household_Mental_Illness      0.34921341
# ACE_Incarcerated_Household_Member 0.99317452
# ACE_Parents_Treated_Violently     0.83591693
# ACE_Score                         0.76848591
# VSI_Score                         0.24496529
# PSS_Score                         0.61844365
# IbsHcClusterMemb                          NA



# between cluster 1 and HC 

map4_4<-subset(map4, map4[,32] == 1 | map4[,32] == 4);dim(map4_4)
traits4<-matrix(NA, nrow=32, ncol=1)
for ( i in 1: 31)
{traits4[i,1]<-summary(lm(map4_4[,i] ~ map4_4[,32]))$coefficients[,4][2]}
row.names(traits4)<-colnames(map4_4)
traits4



# mapping file with 1. IBSHC membership, and IBS membership
map4<-read.delim("C:/Users/swapnajoshi-admin/Documents/swapna/littleIbsMappingIBsHcClust.csv", sep=",", row.names = 1); dim(map4)
map4_ibs<-subset(map4, map4$Group == 2)
map5<-read.delim("C:/Users/swapnajoshi-admin/Documents/swapna/littleIbsMappingIBsClust.csv", sep=",", row.names = 1); dim(map5)
map4_hc<-subset(map4, map4$Group == 1)
map5<-map5[row.names(map4_ibs),]
match(row.names(map5),row.names(map4_ibs))
map4_ibs$IbsClusterMemb<-map5$IbsClusterMemb
map4_hc$IbsClusterMemb<-4
mappingIBS_IBSHC_ClusterMemb<-rbind(map4_ibs, map4_hc)

map4$IBSMembership<-NA
  for ( i in 1:83) {
  if(map4[i,3] == 1) {
    map4[i,34] <- 4
  }
  else if(map4[i,32] != 4) {
    map4[i,34] <- map5[i,34]
  }
}
  for ( i in 1:83) {
    if(map4[i,34] == 4) {
      map4[i,33] <- 4 
    }
    else if(map4[i,32] != 4) {
      map4[i,33] <- map4[i,33]
    }
  }

mappingIBS_IBSHC_ClusterMemb$ACE_Emotional_Abuse<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Sexual_Abuse<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Score<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Physical_Abuse<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Substance_Abuse<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Parental_DivorceSep<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Household_Mental_Illness<-NULL
mappingIBS_IBSHC_ClusterMemb$ACE_Parents_Treated_Violently<-NULL

#write the ibs cluster membership to mapping file

dim(mappingIBS_IBSHC_ClusterMemb)
# [1] 83 26

mappingIBS_IBSHC_ClusterMemb$ibsCluster_1<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,26]==1) {
    mappingIBS_IBSHC_ClusterMemb[i,27]<-1
    }  else (mappingIBS_IBSHC_ClusterMemb[i,27]<-0)
}

mappingIBS_IBSHC_ClusterMemb$ibsCluster_2<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,26]==2) {
    mappingIBS_IBSHC_ClusterMemb[i,28]<-1
    }  else (mappingIBS_IBSHC_ClusterMemb[i,28]<-0)
}

mappingIBS_IBSHC_ClusterMemb$ibsCluster_hc<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,26]==4){
    mappingIBS_IBSHC_ClusterMemb[i,29]<-1
    }  else (mappingIBS_IBSHC_ClusterMemb[i,29]<-0)
}


#write the ibs hc cluster membership to mapping file


mappingIBS_IBSHC_ClusterMemb$ibsHcCluster_1<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,25]==1) {
    mappingIBS_IBSHC_ClusterMemb[i,30]<-1
  }  else (mappingIBS_IBSHC_ClusterMemb[i,30]<-0)
}

mappingIBS_IBSHC_ClusterMemb$ibsHcCluster_2<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,25]==2) {
    mappingIBS_IBSHC_ClusterMemb[i,31]<-1
  }  else (mappingIBS_IBSHC_ClusterMemb[i,31]<-0)
}

mappingIBS_IBSHC_ClusterMemb$ibsHcCluster_3<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,25]==3) {
    mappingIBS_IBSHC_ClusterMemb[i,32]<-1}
  else (mappingIBS_IBSHC_ClusterMemb[i,32]<-0)
}

mappingIBS_IBSHC_ClusterMemb$ibsHcCluster_hc<-NA
for( i in 1: dim(mappingIBS_IBSHC_ClusterMemb)[1])
{
  if(mappingIBS_IBSHC_ClusterMemb[i,25]==4) {
    mappingIBS_IBSHC_ClusterMemb[i,33]<-1
    }  else (mappingIBS_IBSHC_ClusterMemb[i,33]<-0)
}


write.table(mappingIBS_IBSHC_ClusterMemb, file= "mappingIBS_IBSHC_ClusterMemb.csv", sep=',', col.names=NA)



##################### firmicutes bacteroidetes ratio from qiime

fb_ratio <- read.delim("C:/Users/swapnajoshi-admin/Documents/swapna/AlphaBetaDiv04292016/ratios/FirmiBactRatio.txt", sep = "\t", row.names = 1); dim(fb_ratio)

# Group
fb_ratio<-as.matrix(fb_ratio)
fb_ratio<-as.data.frame(fb_ratio[row.names(mappingIBS_IBSHC_ClusterMemb),])
match(row.names(fb_ratio), row.names(mappingIBS_IBSHC_ClusterMemb))

mappingIBS_IBSHC_ClusterMemb$fb_ratio <- fb_ratio[,1]
# summary(lm(mappingIBS_IBSHC_ClusterMemb[,34]~ as.factor(mappingIBS_IBSHC_ClusterMemb[,3])))
# Call:
#   lm(formula = mappingIBS_IBSHC_ClusterMemb[, 34] ~ as.factor(mappingIBS_IBSHC_ClusterMemb[, 
#                                                                                            3]))
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.59543 -0.48184 -0.01313  0.50112  2.05991 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                                    -0.3409     0.1305  -2.611   0.0107 *
#   as.factor(mappingIBS_IBSHC_ClusterMemb[, 3])2  -0.1510     0.1649  -0.916   0.3625  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7268 on 81 degrees of freedom
# Multiple R-squared:  0.01025,	Adjusted R-squared:  -0.001974 
# F-statistic: 0.8385 on 1 and 81 DF,  p-value: 0.3625

boxplot(mappingIBS_IBSHC_ClusterMemb[,34] ~ mappingIBS_IBSHC_ClusterMemb[,3])
p <- ggplot(mappingIBS_IBSHC_ClusterMemb, aes(factor(Group), fb_ratio))
q <- p + geom_boxplot(aes(fill = factor(Group)))
ggsave(q, file="FB_ratio_boxplot_Group.png")


ma.phylum1<-t(ma.phylum)
match(row.names(ma.phylum1), row.names(littleIbsMapping))
littleIbsMapping<-littleIbsMapping[row.names(ma.phylum1),]
ma.phylum1<-as.data.frame(ma.phylum1)
littleIbsMapping$FBRatio<-ma.phylum1$'k__Bacteria; p__Firmicutes'/ma.phylum1$'k__Bacteria; p__Bacteroidetes'
littleIbsMapping$FBRatio<-as.numeric(as.character(littleIbsMapping$FBRatio))
littleIbsMapping<-as.data.frame(littleIbsMapping)

t.test(littleIbsMapping[,31]~ littleIbsMapping[,7])
boxplot(littleIbsMapping[,31]~ littleIbsMapping[,7])

#mean in group HC mean in group IBS 
#       0.9491604         0.7893911  p-value = 0.4082


summary(aov(littleIbsMapping[,31]~ littleIbsMapping[,8]))
littleIbsMapping[,8]<-gsub("U","M",littleIbsMapping[,8])
#p=0.807
boxplot(littleIbsMapping[,31]~ littleIbsMapping[,8])
summary(aov(littleIbsMapping[,31]~ littleIbsMapping[,25]))
#7.17e-07
boxplot(littleIbsMapping[,31]~ littleIbsMapping[,25])


boxplot(littleIbsMapping[,31]~ littleIbsMapping[,7])

t.test(littleIbsMapping[,31]~ littleIbsMapping[,9])
boxplot(littleIbsMapping[,31]~ littleIbsMapping[,9])


t.test(littleIbsMapping[littleIbsMapping$Dx=="IBS",31]~ littleIbsMapping[littleIbsMapping$Dx=="IBS",9])
boxplot(littleIbsMapping[littleIbsMapping$Dx=="IBS",31]~ littleIbsMapping[littleIbsMapping$Dx=="IBS",9])

t.test(littleIbsMapping[littleIbsMapping$Dx=="HC",31]~ littleIbsMapping[littleIbsMapping$Dx=="HC",9])
boxplot(littleIbsMapping[littleIbsMapping$Dx=="HC",31]~ littleIbsMapping[littleIbsMapping$Dx=="HC",9])



##### integration




