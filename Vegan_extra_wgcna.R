

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


