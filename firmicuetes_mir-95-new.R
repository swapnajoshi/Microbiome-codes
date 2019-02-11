#setwd("C:/Users/swapnajoshi/Documents/UCLA_research/IBSmicrobiomeProject/biopsy_scor_microbiome/")
#setwd("C:/Users/swapnajoshi/Documents/UCLA_research/IBSmicrobiomeProject/biopsy_scor_microbiome/")
setwd("C:/Users/swapnajoshi-admin/Documents/swapna/")

install.packages("vegan",  dependencies=T)

library(vegan)

# read OTU table
header1<-scan("otu_table_SCOR_biopsy.txt", nlines = 2, what = character())[8:91]
genus=read.table("otu_table_SCOR_biopsy.txt", sep="\t", skip=2, header=FALSE, row.names = 1)	
names(genus)<-header1

taxa.names = genus$ConsensusLineage

genus1 = as.matrix(genus[,-dim(genus)[2]])
genus1 = scale(genus1, center=F, scale=colSums(genus1))
genus1 = as.data.frame(t(genus1))
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

d.species = otu2taxonomy(genus1,level=7,taxa=taxa.names)
dim(d.species)
d.genus = otu2taxonomy(genus1,level=6,taxa=taxa.names)
dim(d.genus)
d.family = otu2taxonomy(genus1,level=5,taxa=taxa.names)
dim(d.family)
d.order = otu2taxonomy(genus1,level=4,taxa=taxa.names)
dim(d.order)
d.class = otu2taxonomy(genus1,level=3,taxa=taxa.names)
dim(d.class)
d.phylum = otu2taxonomy(genus1,level=2,taxa=taxa.names)
dim(d.phylum)


ma.species = subset(t(d.species), rowMeans(t(d.species)) > 0.01)
dim(ma.species)
 
ma.phylum = subset(t(d.phylum), rowMeans(t(d.phylum)) > 0.01)
dim(ma.phylum)
 

# mapping file


map1 = read.delim("ibs.mapping_new.csv",  sep=',', row.names=1)

map2<-map1[row.names(map1)%in%colnames(ma.species),] ## 2 samples did not have data

map2<-map2[colnames(ma.species),]
match(row.names(map2),colnames(ma.species))
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
# [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
# [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
# [76] 76 77 78 79 80 81 82 83
dim(ma.species)
# [1] 17  83
summary(rowMeans(ma.species))
hist(rowMeans(ma.species))
sorted.means<-rev(sort(rowMeans((ma.species))))
head(sorted.means)

#        k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__ 
                                                                                                          0.31353220 
#         k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__copri 
                                                                                                          0.09533005 
#                           k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__ 
                                                                                                          0.07483813 
#                                          k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__; g__; s__ 
                                                                                                          0.05215287 
#       k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__; s__ 
                                                                                                          0.04639500 
# k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Faecalibacterium; s__prausnitzii 
                                                                                                          0.03567753 


vars<-row.names(ma.phylum)
var.dat<-t(ma.phylum)[row.names(map2),]
par(mar = rep(2, 4))
boxplot(var.dat[,1], map2$Group)

save(d.species, ma.species, map2, genus1, genus, d.phylum, ma.phylum, file="OtuToDataFrame.Rda")

save(ma.species, file="sel.bact.species.matrix.rda")
# Alpha diversitites

num_genera=specnumber(ma.phylum)
chao1=estimateR(ma.phylum*100000000)[2,]
shannon=diversity(ma.phylum,  "shannon")

save(shannon, file="shannon.rda")

# plot alpha diversity

par(mar = rep(2, 4))
boxplot(num_genera  ~	ma.phylum[,1], col=rainbow(3), main="Number of genera")
 boxplot(chao1~  ma.phylum,  col=rainbow(3),	 main="Chao1  index")
boxplot(shannon~  ma.phylum[,1],	  col=rainbow(3),  main="Shannon  evenness")


# summary(aov(shannon  ~  enterotype))


# Comparison of bacterial communities
# 1. which species are different between ibs and hc,
# 2. diffrence between bowel habbit subtypes


# heatmap

#########################################################
setwd("C:/Users/swapnajoshi/Documents/UCLA_research/IBSmicrobiomeProject/biopsy_scor_microbiome/")
load("OtuToDataFrame.Rda")
#WGCNA
library(WGCNA)
library(flashClust)
#Heatmap
library(dplyr)
library(NMF)
library(RColorBrewer)

# trait data 
traitData = map2 #from above

dim(ma.species)
ma.species[1:2,1:2]
datExpr0<-t(ma.species)
#datExpr0<-subset(ma.phylum, map2$Dx=="IBS")
#map2.ibs<-subset(map2, map2$Dx=="IBS")
# row.names(datExpr0)<-substr(row.names(datExpr0),8,13)


gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

sampleTree = flashClust(dist(datExpr0), method = "complete")
match(row.names(map2),row.names(datExpr0)) ## ok
labels1=map2$Dx
dat.ordered<-datExpr0[sampleTree$order,]

# otuData.ibs<-subset(datExpr0, row.names(datExpr0)%in% [traitData grep("IBS",colnames(otuData)),])
# otuData.hc<-datExpr0[grep("HC",colnames(otuData)),]

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers_new", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
abline(h = 0.55, col = "red")

clust1 = cutreeStatic(sampleTree, cutHeight = 0.55, minSize = 1)
table(clust1)
# clust1

#clust1
# 1  2  3  4 
#44 20 12  7 

colStates<-factor(row.names(dat.ordered),levels=sampleTree$labels[sampleTree$order])
Clust2<-clust1[sampleTree$order]
lab<-map2[sampleTree$labels,]$Dx
df2<-as.data.frame(cbind(Clust2,as.character(lab)))
head(df2)
set.seed(123)
rainbow(n=5)
#[1] "#FF0000FF" "#00FF00FF" "#0000FFFF" "#0066FFFF" "#CC00FFFF"

df2$color1<-NA
for ( i in 1: 83)
{if (df2[i,1]=="1")
	{df2[i,3]<-"#FF0000FF"}
else if (df2[i,1]=="2")
	{df2[i,3]<-"#00FF00FF"}
else if (df2[i,1]=="3")
	{df2[i,3]<-"#0000FFFF"}
}

df2$colorDx<-NA
for ( i in 1: 83)
{if (df2[i,2]=="IBS")
	{df2[i,4]<-"#0066FFFF"}
else if (df2[i,2]=="HC")
	{df2[i,4]<-"#CC00FFFF"}
}

plot(sampleTree, main = "Sample clustering to detect outliers_new", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, labels=labels1)
abline(h = 0.55, col = "red")
colored_bars(df2[,c(3,4)], sampleTree)



#############################################
  
memb1<-as.matrix(cbind(rownames(datExpr0),clust1))
memb1<-as.data.frame(memb1)
memb1$labels<-labels1
table(memb1[,2],memb1[,3])
#  HC IBS
#  1 19  25
#  2  6  14
#  3  5   7
#  4  1   6
exp.clu1<-t(datExpr0[row.names(datExpr0)%in%memb1[memb1[,2]=="1" & memb1[,3]!="HC",1],])
exp.clu2<-t(datExpr0[row.names(datExpr0)%in%memb1[memb1[,2]=="2" & memb1[,3]!="HC",1],])
exp.clu3<-t(datExpr0[row.names(datExpr0)%in%memb1[memb1[,2]=="3" & memb1[,3]!="HC",1],])
exp.clu4<-t(datExpr0[row.names(datExpr0)%in%memb1[memb1[,2]=="4" & memb1[,3]!="HC",1],])
exp.hc<-t(datExpr0[row.names(datExpr0)%in%memb1[memb1[,3]=="HC",1],])


save(exp.clu1, file="exp.clu1.rda")
save(exp.clu2, file="exp.clu2.rda")
save(exp.clu3, file="exp.clu3.rda")
save(exp.clu4, file="exp.clu4.rda")
save(exp.hc, file="exp.hc.rda")

########################################################################

row.names(memb1)<-memb1[,1]
map1<-merge(map2, memb1, by="row.names")
row.names(map1)<-map1[,1]
map1$V1<-NULL
map1$labels<-NULL
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
if(map1[i,79]==3){map1[i,82]<-1}
else if(map1[i,79]!=3) {map1[i,82]<-0}
}

map1$cluster4<-NA
for( i in 1: dim(df2)[1])
{
if(map1[i,79]==4){map1[i,83]<-1}
else if(map1[i,79]!=4) {map1[i,83]<-0}
}

map1$HC<-NA
for( i in 1: dim(df2)[1])
{
if(map1[i,79]=="HC"){map1[i,84]<-1}
else if(map1[i,79]!="HC") {map1[i,84]<-0}
}

map2<-map1[,c(1:77,79:84,78)]
colnames(map2)[78]<-"cluster"


write.table(map2, file="bigIbsMapping.csv",sep=",", col.names=NA)
colnames(map2)
littleIbsMapping<-map2[,c(1,2,3,4,5,6,8,9,10,30:42,76:83)]
write.table(littleIbsMapping, file="littleIbsMapping.csv",sep=",", col.names=NA)


#####################################################################
col.j<-jet.colors(75)
col.mat<-as.matrix(df2[,c(3,4)])
png("heatmap.ibs.hc.complete.png", bg="white", res=300, width=3000, height=3000)
hv2<-heatmap.plus(as.matrix(t(dat.ordered)),na.rm=TRUE,scale="none",  
               #RowSideColor=col.mat,
                ColSideColors=col.mat,
              col=col.j,
               #		zlim=c(0,1),
               Colv=NA,Rowv=T,
               cexRow=1,cexCol=1,
               #		hclustfun = dist1,
               # keep.dendro = NA,
               main = paste("IBS_HC", dim(t(dat.ordered))[1],"Probes; ",dim(t(dat.ordered))[2],"Samples"),
               labRow=row.names(t(dat.ordered))
)
dev.off()


####################################################################################

####################################################################################
littleIbsMapping<-read.delim("littleIbsMapping.csv", sep=',', row.names=1)
head(littleIbsMapping)
load("exp.clu1.rda")
load("exp.clu2.rda")
load("exp.clu3.rda")
load("exp.clu4.rda")
load("exp.hc.rda")
datExpr1<-datExpr0[row.names(littleIbsMapping),]
datExpr1<-datExpr1[row.names(littleIbsMapping[littleIbsMapping$Dx!="HC",]),]
littleIbsMapping1<-littleIbsMapping[littleIbsMapping$Dx!="HC",]
match(row.names(datExpr1), row.names(littleIbsMapping1))
littleIbsMapping1<-as.data.frame(littleIbsMapping1)
littleIbsMapping1$BarcodeSequence <-as.numeric(as.factor(littleIbsMapping1$BarcodeSequence))
littleIbsMapping1<-littleIbsMapping1[,-c(1,2,3,4,5,6)]
littleIbsMapping1$Dx<-gsub("IBS", 2, littleIbsMapping1$Dx)
littleIbsMapping1$Dx<-gsub("HC", 1, littleIbsMapping1$Dx)
littleIbsMapping1$Bowel.Habit<-gsub("C", 1, littleIbsMapping1$Bowel.Habit)

oneway.test(as.numeric(littleIbsMapping1[,19])~littleIbsMapping1[,10])$p.value
littleIbsMapping


table(littleIbsMapping1$cluster)
mic.df<-as.data.frame(t(cbind(exp.clu1,exp.clu2)))
mic.df$memb<-c(rep("1",25),rep("2",14))
ttest.dis.cont<-t(sapply(mic.df[-18], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
ttest.dis.cont$FC<-NA
ttest.dis.cont.sig<-subset(ttest.dis.cont, bh<0.05)
save(ttest.dis.cont, file="ttest.clu1.clu2.rda")
write.table(ttest.dis.cont, file="ttest.clu1.clu2.csv", sep=",", col.names=NA)

mic.df<-as.data.frame(t(cbind(exp.clu1,exp.clu3)))
mic.df$memb<-c(rep("1",25),rep("2",7))
ttest.dis.cont<-t(sapply(mic.df[-18], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
ttest.dis.cont$FC<-NA
ttest.dis.cont.sig<-subset(ttest.dis.cont, bh<0.05)
save(ttest.dis.cont, file="ttest.clu1.clu3.rda")
write.table(ttest.dis.cont, file="ttest.clu1.clu3.csv", sep=",", col.names=NA)


mic.df<-as.data.frame(t(cbind(exp.clu1,exp.clu4)))
mic.df$memb<-c(rep("1",25),rep("2",6))
ttest.dis.cont<-t(sapply(mic.df[-18], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
ttest.dis.cont$FC<-NA
ttest.dis.cont.sig<-subset(ttest.dis.cont, bh<0.05)
save(ttest.dis.cont, file="ttest.clu1.clu4.rda")
write.table(ttest.dis.cont, file="ttest.clu1.clu4.csv", sep=",", col.names=NA)

mic.df<-as.data.frame(t(cbind(exp.clu2,exp.clu3)))
mic.df$memb<-c(rep("1",14),rep("2",7))
ttest.dis.cont<-t(sapply(mic.df[-18], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
ttest.dis.cont$FC<-NA
ttest.dis.cont.sig<-subset(ttest.dis.cont, bh<0.05)
save(ttest.dis.cont, file="ttest.clu2.clu3.rda")
write.table(ttest.dis.cont, file="ttest.clu2.clu3.csv", sep=",", col.names=NA)


mic.df<-as.data.frame(t(cbind(exp.clu2,exp.clu4)))
mic.df$memb<-c(rep("1",14),rep("2",6))
ttest.dis.cont<-t(sapply(mic.df[-18], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
ttest.dis.cont$FC<-NA
ttest.dis.cont.sig<-subset(ttest.dis.cont, bh<0.05)
save(ttest.dis.cont, file="ttest.clu2.clu4.rda")
write.table(ttest.dis.cont, file="ttest.clu2.clu4.csv", sep=",", col.names=NA)


mic.df<-as.data.frame(t(cbind(exp.clu3,exp.clu4)))
mic.df$memb<-c(rep("1",7),rep("2",6))
ttest.dis.cont<-t(sapply(mic.df[-18], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
ttest.dis.cont<-as.data.frame(ttest.dis.cont)
ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
ttest.dis.cont.sig<-subset(ttest.dis.cont, bh<0.05)
save(ttest.dis.cont, file="ttest.clu3.clu4.rda")
write.table(ttest.dis.cont, file="ttest.clu3.clu4.csv", sep=",", col.names=NA)



littleIbsMapping1$Gender<-gsub("M","1",littleIbsMapping1$Gender)
littleIbsMapping1$Gender<-gsub("F","2",littleIbsMapping1$Gender)
map4<- subset(littleIbsMapping1,  littleIbsMapping1$cluster=="1"| littleIbsMapping1$cluster=="2")
map5<-map4[,-c(1,2)]
map5[,17]<-as.numeric(as.character(map5[,17]))
save(map5, file="map.numeric.clu1.2.rda")
p.traits1<-matrix(NA, nrow=16, ncol=3)
for ( i in 1: 16)
{p.traits1[i,1]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["p.value"])
p.traits1[i,c(2,3)]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["estimate"])
}
row.names(p.traits1)<-colnames(map5)[1:16]
 

map4<- subset(littleIbsMapping1,  littleIbsMapping1$cluster=="1"| littleIbsMapping1$cluster=="3")
map5<-map4[,-c(1,2)]
save(map5, file="map.numeric.clu1.3.rda")
p.traits2<-matrix(NA, nrow=16, ncol=3)
for ( i in 1: 16)
{p.traits2[i,1]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["p.value"])
p.traits2[i,c(2,3)]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["estimate"])
}
row.names(p.traits2)<-colnames(map5)[1:16]


map4<- subset(littleIbsMapping1,  littleIbsMapping1$cluster=="1"| littleIbsMapping1$cluster=="4")
map5<-map4[,-c(1,2)]
save(map5, file="map.numeric.clu1.4.rda")
p.traits3<-matrix(NA, nrow=16, ncol=3)
for ( i in 1: 16)
{for ( i in 1: 16)
{p.traits3[i,1]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["p.value"])
p.traits3[i,c(2,3)]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["estimate"])
}
row.names(p.traits3)<-colnames(map5)[1:16]

map4<- subset(littleIbsMapping1,  littleIbsMapping1$cluster=="2"| littleIbsMapping1$cluster=="3")
map5<-map4[,-c(1,2)]
save(map5, file="map.numeric.clu2.3.rda")
p.traits4<-matrix(NA, nrow=16, ncol=1)
for ( i in 1: 16)
{p.traits4[i,1]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["p.value"])}
row.names(p.traits4)<-colnames(map5)[1:16]


map4<- subset(littleIbsMapping1,  littleIbsMapping1$cluster=="2"| littleIbsMapping1$cluster=="4")
map5<-map4[,-c(1,2)]
save(map5, file="map.numeric.clu2.4.rda")
p.traits5<-matrix(NA, nrow=16, ncol=1)
for ( i in 1: 16)
{p.traits5[i,1]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["p.value"])}
row.names(p.traits5)<-colnames(map5)[1:16]


map4<- subset(littleIbsMapping1,  littleIbsMapping1$cluster=="3"| littleIbsMapping1$cluster=="4")
map5<-map4[,-c(1,2)]
save(map5, file="map.numeric.clu3.4.rda")
p.traits6<-matrix(NA, nrow=16, ncol=1)
for ( i in 1: 16)
{p.traits6[i,1]<-unlist(t.test(as.numeric(as.character(map5[,i]))~map5$cluster, na.rm=TRUE)["p.value"])}
row.names(p.traits6)<-colnames(map5)[1:16]


p.traits<-cbind(p.traits1,p.traits2,p.traits3,p.traits4,p.traits5,p.traits6)


write.table(p.traits, file="t-test-pheno-groupsibs.csv", col.names=NA, sep=",")
                         
match(row.names(datExpr0), row.names(littleIbsMapping))
littleIbsMapping<-littleIbsMapping[row.names(datExpr0),]   
 p.Dx.gender<-matrix(NA, nrow=17, ncol=1)
for ( i in 1:17)                   
{p.Dx.gender[i,1]<-summary(lm(datExpr0[,i]~ littleIbsMapping[,7] * littleIbsMapping[,9]))$coefficients[,4][2]
}

row.names(p.Dx.gender)<-colnames(datExpr0)

write.table(p.Dx.gender, file="p.Dx.gender.csv", col.names=NA)

sig.gen.int<-cbind(datExpr0[,c(4,9)],cbind(as.character(littleIbsMapping[,7]) ,as.character(littleIbsMapping[,9])))
colnames(sig.gen.int)<-c("Bacteroidetes_distasonis","Firmicutes_Lachnospiraceae","Dx","Gender")
sig.gen.int<-as.data.frame(sig.gen.int)
sig.gen.int$int1<-interaction(sig.gen.int$Dx,sig.gen.int$Gender)
ggplot(aes(y = as.numeric(as.character(Bacteroidetes_distasonis)), x = int1), data = sig.gen.int) + geom_boxplot()
ggplot(aes(y = as.numeric(as.character(Firmicutes_Lachnospiraceae)), x = int1), data = sig.gen.int) + geom_boxplot()


a1<-subset(sig.gen.int, sig.gen.int$int1=="HC.F"|sig.gen.int$int1=="IBS.F")
a1[,1]<-as.numeric(as.character(a1[,1]))
a1[,2]<-as.numeric(as.character(a1[,2]))
a1<-as.data.frame(a1)
t.test(a1[,1]~as.character(a1[,5]))$p.value
#p=0.2
t.test(a1[,2]~as.character(a1[,5]))$p.value
#p=0.09

a1<-subset(sig.gen.int, sig.gen.int$int1=="HC.F"|sig.gen.int$int1=="HC.M")
a1[,1]<-as.numeric(as.character(a1[,1]))
a1[,2]<-as.numeric(as.character(a1[,2]))
a1<-as.data.frame(a1)
t.test(a1[,1]~as.character(a1[,5]))$p.value
#p=0.22
t.test(a1[,2]~as.character(a1[,5]))$p.value
#p=0.2

a1<-subset(sig.gen.int, sig.gen.int$int1=="HC.F"|sig.gen.int$int1=="IBS.M")
a1[,1]<-as.numeric(as.character(a1[,1]))
a1[,2]<-as.numeric(as.character(a1[,2]))
a1<-as.data.frame(a1)
t.test(a1[,1]~as.character(a1[,5]))$p.value
#p=0.48
t.test(a1[,2]~as.character(a1[,5]))$p.value
#p=0.09

a1<-subset(sig.gen.int, sig.gen.int$int1=="IBS.F"|sig.gen.int$int1=="HC.M")
a1[,1]<-as.numeric(as.character(a1[,1]))
a1[,2]<-as.numeric(as.character(a1[,2]))
a1<-as.data.frame(a1)
t.test(a1[,1]~as.character(a1[,5]))$p.value
#p=0.93
t.test(a1[,2]~as.character(a1[,5]))$p.value
#p=0.61


a1<-subset(sig.gen.int, sig.gen.int$int1=="IBS.F"|sig.gen.int$int1=="IBS.M")
a1[,1]<-as.numeric(as.character(a1[,1]))
a1[,2]<-as.numeric(as.character(a1[,2]))
a1<-as.data.frame(a1)
t.test(a1[,1]~as.character(a1[,5]))$p.value
#p=0.28
t.test(a1[,2]~as.character(a1[,5]))$p.value
#p=0.79


a1<-subset(sig.gen.int, sig.gen.int$int1=="HC.M"|sig.gen.int$int1=="IBS.M")
a1[,1]<-as.numeric(as.character(a1[,1]))
a1[,2]<-as.numeric(as.character(a1[,2]))
a1<-as.data.frame(a1)
t.test(a1[,1]~as.character(a1[,5]))$p.value
#p=0.3
t.test(a1[,2]~as.character(a1[,5]))$p.value
#p=0.5


##################### firmicutes bacteroidetes ratio at phylum level

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

# import miRNA data:
  # import microbiome data
# make the matrices with same samples and same order

load("C:/Users/swapnajoshi/Documents/UCLA_research/miRNAtissueSerumIbs/nanostring_data/raw_data/all_rccs/TissueRCC/norm.data4.rda")
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
# [38] 38 39 40 41 42 43
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
ctab <- cor(data1, method="spearman")

#library(ltm)
#ctab.p <- rcor.test(data1, method="spearman")

corMat<-round(ctab, 2)
# sub1<-function(x)(subset(x, x[i]>=0.5))
highCor<-matrix(NA, nrow=507, ncol=507)
row.names(highCor)<-row.names(corMat)
colnames(highCor)<-colnames(corMat)
for ( i in 1:dim(corMat)[1])
  for ( j in 1:dim(corMat)[2])
  {
    if(corMat[i,j]>0.5)
    {highCor[i,j]<- 1}
    else if(corMat[i,j]< -0.5)
    {highCor[i,j]<- 1}
    else (highCor[i,j]<- 0)
  }

mirnaMicroCor<-highCor[1:490,491:507]
mirnaMicroCor<-as.data.frame(mirnaMicroCor)

mirnaMicroCor$sum1<-apply(mirnaMicroCor,1,sum)

mirnaMicroCor1<-mirnaMicroCor[mirnaMicroCor$sum1>0,]
data1<-as.data.frame(data1)
dim(mirnaMicroCor1)

write.table(mirnaMicroCor1, file="mirnaMicroCor1.csv", sep=",", col.names=NA)


cor.microb<-t(microbMat2[c( "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Porphyromonadaceae; g__Parabacteroides; s__distasonis","k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__; s__" ,"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Coprococcus; s__" ),])
cor.mirna<-t(mirna1[c( "hsa-miR-27b-3p|0"   ,   "hsa-miR-95|0"    ,      "hsa-miR-1|0"     ,      "hsa-miR-145-5p|0"  ,   
 "hsa-miR-196a-5p|0"  ,   "hsa-miR-99b-5p|0"    ,  "hsa-miR-125b-5p|0.009", "hsa-miR-331-3p|0"     ,
 "hsa-miR-139-5p|0"   ,   "hsa-miR-125a-5p|0"  ),])

data2<-as.data.frame(cbind(cor.microb,cor.mirna))
map.mirna<-littleIbsMapping[substr(row.names(littleIbsMapping),8,13)%in%row.names(data2),]
dim(map.mirna)
row.names(map.mirna)<-substr(row.names(map.mirna),8,13)
map.mirna<-map.mirna[row.names(data2),]
match(row.names(data2),row.names(map.mirna))
data3<-cbind(data2, map.mirna[,c(7,8,9,25)])
colnames(data3)<-gsub(" ","",colnames(data3))
colnames(data3)<-gsub(";","",colnames(data3))
colnames(data3)<-gsub("-","_",colnames(data3))
colnames(data3)<-gsub("\\|","",colnames(data3))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_27b_3p0))
P+geom_point(aes(color=Dx))
ggsave(firmicutes_coprococcus_
cor.test(data3$k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,data3$hsa_miR_27b_3p0, method='spearman')
#p.value  5.197e-05
P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_27b_3p0))
P+geom_point(aes(color=Bowel.Habit))
cor.test(data3$k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,data3$hsa_miR_27b_3p0,  method='spearman')

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_27b_3p0))
P+geom_point(aes(color=as.factor(cluster)))



P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_950))
P+geom_point(aes(color=Dx, size=4))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_950))
P+geom_point(aes(color=Bowel.Habit))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_950))
P+geom_point(aes(color=as.factor(cluster), size=4))


P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_27b_3p0))
P+geom_point(aes(color=Dx, size=4))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_27b_3p0))
P+geom_point(aes(color=Bowel.Habit))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_27b_3p0))
P+geom_point(aes(color=as.factor(cluster), size=4))




P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_99b_5p0))
P+geom_point(aes(color=Dx, size=4))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_99b_5p0))
P+geom_point(aes(color=Bowel.Habit))

P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_99b_5p0))
P+geom_point(aes(color=as.factor(cluster), size=4))





P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_950,color=as.factor(Dx)))
P+geom_point(size=4)+
scale_colour_hue(l=50) +
geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=TRUE)



P<-ggplot(data3, aes(k__Bacteriap__Firmicutesc__Clostridiao__Clostridialesf__Lachnospiraceaeg__Coprococcuss__,hsa_miR_950,color=as.factor(cluster)))
P+geom_point(size=4)+
#scale_colour_hue(l=50) +
geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE)








data2<-rbind(mirna[c("hsa-let-7g-5p|0","hsa-miR-34c-5p|0","hsa-miR-451a|0","hsa-miR-25-3p|0","hsa-miR-374b-5p|0"),], microbMat2["k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Porphyromonadaceae; g__Parabacteroides; s__distasonis",])
data2<-t(data2)
data2<-rbind( mirna["k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Coprococcus; s__",], microbMat2["hsa-miR-99b-5p|0",])
plot(data2$"k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Coprococcus; s__",data2$"hsa-miR-99b-5p|0")

colnames(data2)[6]<-"distenosis"
data2<-as.data.frame(data2)
plot(data2$"hsa-let-7g-5p|0",data2$distenosis)
plot(data2$"hsa-miR-34c-5p|0",data2$proteobacteria)
plot(data2$"hsa-miR-451a|0",data2$proteobacteria)










# predominant ibs vs normal
exp.hc<-as.data.frame(exp.hc)
exp.hc$memb<-c(rep("Healthy control",31))
mic.df<-rbind(ibs.pre,exp.hc)
ttest.dis.cont<-t(sapply(mic.df[-8], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
                  ttest.dis.cont<-as.data.frame(ttest.dis.cont)
                  ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
                  # ttest.dis.cont$logFC<-ttest.dis.cont$'estimate.mean in group Relapse'-ttest.dis.cont$'estimate.mean in group No Relapse'
                  # ttest.dis.cont$abslogFC<-abs(ttest.dis.cont$logFC)
                  
                  
                  # ttest.dis.cont$absFC<-2^ttest.dis.cont$abslogFC
                  ttest.dis.cont$FC<-NA
                  
                  for (i in 1: 7) {
                    if(ttest.dis.cont[i,2]/ttest.dis.cont[i,1]>1)
                    {ttest.dis.cont[i,8]<-ttest.dis.cont[i,2]/ttest.dis.cont[i,1]}
                    else if (ttest.dis.cont[i,2]/ttest.dis.cont[i,1]<=1)
                    {ttest.dis.cont[i,8]<--1*(ttest.dis.cont[i,1]/ttest.dis.cont[i,2])}
                  }
                  
                  save(ttest.dis.cont, file="ttest.ibs.vs.controls.rda")
                  write.table(ttest.dis.cont, file="ttest.ibs.vs.controls.csv", sep=",", col.names=NA)
                  
                  
                  
 # normal-like ibs vs control
                  exp.hc<-as.data.frame(exp.hc)
                  exp.hc$memb<-c(rep("Healthy control",31))
                  mic.df<-rbind(ibs.norLike,exp.hc)
                  ttest.dis.cont<-t(sapply(mic.df[-8], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
                  ttest.dis.cont<-as.data.frame(ttest.dis.cont)
                  ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
                  # ttest.dis.cont$logFC<-ttest.dis.cont$'estimate.mean in group Relapse'-ttest.dis.cont$'estimate.mean in group No Relapse'
                  # ttest.dis.cont$abslogFC<-abs(ttest.dis.cont$logFC)
                  
                  
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


# Genus level differences
                  d.genus= otu2taxonomy(genus1,level=5,taxa=taxa.names)  
                  genus.ibs<-d.genus[row.names(d.genus)%in%memb1[,1],]
                  genus.hc<-d.genus[!row.names(d.genus)%in%memb1[,1],]
                 match( row.names(genus.ibs),memb1[,1])
                  genus.df<-cbind(genus.ibs, memb1[,2])
                  ibs.pre<-genus.ibs[row.names(genus.ibs)%in%memb1[memb1[,2]!=2,][,1],]
                  ibs.norLike<-genus.ibs[row.names(genus.ibs)%in%memb1[memb1[,2]==2,][,1],]
                  exp.hc<-d.genus[!row.names(d.genus)%in%memb1[,1],]
                  ibs.pre<-as.data.frame(ibs.pre)
                  ibs.pre$memb<-c(rep("ibs_predominant",34))
                  ibs.norLike<-as.data.frame(ibs.norLike)
                  ibs.norLike$memb<-c(rep("normal_like",18))
                  
                  mic.df<-rbind(ibs.pre)
                  ttest.dis.cont<-t(sapply(mic.df[-189], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
                                    ttest.dis.cont<-as.data.frame(ttest.dis.cont)
                                    ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
                                    # ttest.dis.cont$logFC<-ttest.dis.cont$'estimate.mean in group Relapse'-ttest.dis.cont$'estimate.mean in group No Relapse'
                                    # ttest.dis.cont$abslogFC<-abs(ttest.dis.cont$logFC)
                                    
                                    ttest.dis.cont.sig<-subset(ttest.dis.cont,ttest.dis.cont$bh<=0.05 & ttest.dis.cont$p.value!= "NaN")
                                    # ttest.dis.cont$absFC<-2^ttest.dis.cont$abslogFC
                                    ttest.dis.cont.sig$family<-unlist(strsplit(row.names(ttest.dis.cont.sig),";"))[grep("f__",unlist(strsplit(row.names(ttest.dis.cont.sig),";")))]
                                    ttest.dis.cont.sig$FC<-NA
                                    
                                    for (i in 1: 5) {
                                      if(ttest.dis.cont.sig[i,2]/ttest.dis.cont.sig[i,1]>1)
                                      {ttest.dis.cont.sig[i,9]<-ttest.dis.cont.sig[i,2]/ttest.dis.cont.sig[i,1]}
                                      else if (ttest.dis.cont.sig[i,2]/ttest.dis.cont.sig[i,1]<=1)
                                      {ttest.dis.cont.sig[i,9]<--1*(ttest.dis.cont.sig[i,1]/ttest.dis.cont.sig[i,2])}
                                    }
                                    
                                    save(ttest.dis.cont.sig, file="ttest.dis.cont.sig.family_level_pre_norlike.rda")
                                    write.table(ttest.dis.cont.sig, file="ttest.dis.cont.sig.family_level_pre_norlike.rda.csv", sep=",", col.names=NA)
                                    
                                    
                                    # predominant ibs vs normal
                                    exp.hc<-as.data.frame(exp.hc)
                                    exp.hc$memb<-c(rep("Healthy control",31))
                                    mic.df<-rbind(ibs.pre,exp.hc)
                                    ttest.dis.cont<-t(sapply(mic.df[-8], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
                                    ttest.dis.cont<-as.data.frame(ttest.dis.cont)
                                    ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
                                    # ttest.dis.cont$logFC<-ttest.dis.cont$'estimate.mean in group Relapse'-ttest.dis.cont$'estimate.mean in group No Relapse'
                                    # ttest.dis.cont$abslogFC<-abs(ttest.dis.cont$logFC)
                                    
                                    
                                    # ttest.dis.cont$absFC<-2^ttest.dis.cont$abslogFC
                                    ttest.dis.cont$FC<-NA
                                    
                                    for (i in 1: 7) {
                                      if(ttest.dis.cont[i,2]/ttest.dis.cont[i,1]>1)
                                      {ttest.dis.cont[i,8]<-ttest.dis.cont[i,2]/ttest.dis.cont[i,1]}
                                      else if (ttest.dis.cont[i,2]/ttest.dis.cont[i,1]<=1)
                                      {ttest.dis.cont[i,8]<--1*(ttest.dis.cont[i,1]/ttest.dis.cont[i,2])}
                                    }
                                    
                                    save(ttest.dis.cont, file="ttest.ibs.vs.controls.rda")
                                    write.table(ttest.dis.cont, file="ttest.ibs.vs.controls.csv", sep=",", col.names=NA)
                                    
                                    
                                    
                                    # normal-like ibs vs control
                                    exp.hc<-as.data.frame(exp.hc)
                                    exp.hc$memb<-c(rep("Healthy control",31))
                                    mic.df<-rbind(ibs.norLike,exp.hc)
                                    ttest.dis.cont<-t(sapply(mic.df[-8], function(x) 	unlist(t.test(x~mic.df$memb, na.rm=TRUE)[c("estimate","p.value","statistic","conf.int")])))
                                    ttest.dis.cont<-as.data.frame(ttest.dis.cont)
                                    ttest.dis.cont$bh<-p.adjust(ttest.dis.cont$"p.value", method='fdr')
                                    # ttest.dis.cont$logFC<-ttest.dis.cont$'estimate.mean in group Relapse'-ttest.dis.cont$'estimate.mean in group No Relapse'
                                    # ttest.dis.cont$abslogFC<-abs(ttest.dis.cont$logFC)
                                    
                                    
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
                                    

                                    
       # Pi chart family
                                    
                                    png("pie.chart.family.png", res=300, width=3000, height=1500)
                                    
                                    par(mfrow=  c(1,3), mar=  c(1, 2, 2, 8))
                                    slices.ibs.pre<-round(colMeans(ibs.pre.sig), digits=0)
                                    
                                    lbls <- names(slices.ibs.pre)
                                    pct <- round(slices.ibs.pre/sum(slices.ibs.pre)*100)
                                    # lbls <- paste(lbls, pct) # add percents to labels
                                    # lbls <- paste(lbls,"%",sep="") # ad % to labels
                                    pie(pct,labels = NA, col=rainbow(length(lbls)),radius = 1.0,
                                        main="Pie Chart IBS") 
                                    
                                    
                                    
                                    slices.ibs.norLike<-round(colMeans(ibs.norLike.sig), digits=0)
                                    
                                    # slices <- c(10, 12, 4, 16, 8)
                                    
                                    
                                    lbls <- names(slices.ibs.norLike)
                                    pct <- round(slices.ibs.norLike/sum(ibs.norLike.sig)*100)
                                    # lbls <- paste(lbls, pct) # add percents to labels
                                    # lbls <- paste(lbls,"%",sep="") # ad % to labels
                                    pie(pct,labels = NA, col=rainbow(length(lbls)),radius = 1.0,
                                        main="Pie Chart Normal-like IBS") 
                                    
                                    
                                    slices.hc<-round(colMeans(exp.hc.sig), digits=0)
                                    
                                    # slices <- c(10, 12, 4, 16, 8)
                                    
                                    
                                    lbls <- names(slices.hc)
                                    pct <- round(slices.hc/sum(slices.hc)*100)
                                    # lbls <- paste(lbls, pct) # add percents to labels
                                    # lbls <- paste(lbls,"%",sep="") # ad % to labels
                                    pie(pct,labels = NA, col=rainbow(length(lbls)),radius = 1.0,
                                        main="Pie Chart healthy controls") 
                                    
                                    legend("topleft", c("p__Bacteroidetes;f__Bacteroidaceae",  "p__Bacteroidetes;f__Porphyromonadaceae" ,  "p__Firmicutes; o__Clostridiales; f__" , 
                                                         "p__Firmicutes; f__Lachnospiraceae"  , "p__Firmicutes; f__Ruminococcaceae",  "k__Bacteria; p__Tenericutes"),  cex=1, fill=rainbow(length(lbls)))
                                    dev.off()
                                    
                                    
                                    phylum.ibs.pre<-ma.phylum[row.names(ma.phylum)%in%colnames(ibs.pre),]
                                    phylum.ibs.nlike<-ma.phylum[row.names(ma.phylum)%in%colnames(ibs.pre),]
                                    map2.ibs.pre<-map2.ibs[row.names(map2.ibs)%in%colnames(ibs.pre),]

                                   match(row.names(map2.ibs.pre),row.names(phylum.ibs.pre))
                                   df1<-cbind(phylum.ibs.pre[,-c(7,16,19,22)],map2.ibs.pre[,c(8,9,10,11,21:42)])
                                   library(ggplot2)
                                   library(reshape2)
                                   qplot(x=Var1, y=Var2, data=melt(cor(df1,use = 'pairwise.complete.obs')), fill=value, geom="tile")
                                   
                                   df2<-df1[-8,]
                                           cor.test(df1$ACE_Sexual_Abuse,df1$'k__Bacteria; p__Proteobacteria', use = 'pairwise.complete.obs')
                                    
#                                    ACE_Physical_Abuse    k__Bacteria; p__Elusimicrobia 0.4651750
#                                    HAD_Depression     k__Bacteria; p__Fusobacteria 0.4039963
#                                    ACE_Physical_Abuse    k__Bacteria; p__Lentisphaerae 0.4250854
#                                    ACE_Physical_Abuse   k__Bacteria; p__Proteobacteria 0.6138615
#                                    ACE_Physical_Abuse     k__Bacteria; p__Spirochaetes 0.4324277
#                                    
                                   plot(df2$HAD_Depression,df2$'k__Bacteria; p__Fusobacteria')
                                   
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
row.names(microbMat)<-substr(row.names(microbMat),8,12)
microbMat1<-t(microbMat)
microbMat2<-microbMat1[,colnames(microbMat1)%in%colnames(mirna)]
mirna1<-mirna[,colnames(mirna)%in%colnames(microbMat2)]
microbMat2<-microbMat2[,colnames(mirna1)]
match(colnames(microbMat2),colnames(mirna1))
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
# [38] 38 39 40 41 42 43
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
ctab <- cor(data1)
corMat<-round(ctab, 2)
# sub1<-function(x)(subset(x, x[i]>=0.5))
highCor<-matrix(NA, nrow=545, ncol=545)
row.names(highCor)<-row.names(corMat)
colnames(highCor)<-colnames(corMat)
for ( i in 1:dim(corMat)[1])
  for ( j in 1:dim(corMat)[2])
  {
    if(corMat[i,j]>0.7)
    {highCor[i,j]<- 1}
    else if(corMat[i,j]< -0.7)
    {highCor[i,j]<- 1}
    else (highCor[i,j]<- 0)
  }

mirnaMicroCor<-highCor[,491:545]
mirnaMicroCor<-as.data.frame(mirnaMicroCor)

mirnaMicroCor$sum1<-apply(mirnaMicroCor,1,sum)

mirnaMicroCor1<-mirnaMicroCor[mirnaMicroCor$sum1>0,]
data1<-as.data.frame(data1)

plot(data1$'hsa-miR-34c-5p|0', data1$'k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Sutterella')
