#1. PCoA:  IBS only
#2. PCoA: without the IBS-U- is this referring to the IBS only ? I don’t think you need to remove IBS-U unless you are specifically evaluating bowel habit variable
#3. Association of subtypes with Dx, BH, Sex, and other clinical charachteristics (IBS + HC; IBS only)
#4. FDR for dietary correlations--- done
#5. Correlation of microbial subtypes (IBS+HC; IBS only) with diet
#6. Correlation of microbiome with the new diet variable (Mediterranean, standard and modified American, Vegetarian/Vegan, gluten free)

#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")


# Batch effects correction didnot work well : proceed without correcting
##########################################################################
library(ggplot2)

ps2 <- prune_samples(sample_names(ps) != "A5807", ps)
ps2.ibs <- prune_samples(sample_data(ps2)$Dx == "IBS", ps2) 
ps2.hc <- prune_samples(sample_data(ps2)$Dx == "HC", ps2) 

alpha.diversity <- estimate_richness (ps2.ibs, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps2.ibs), alpha.diversity)

dx.aov <- aov(Chao1 ~ Sex + SequencingBatch, data1) #0.68
dx.aov <- aov(Simpson ~ Sex + SequencingBatch, data1) #0.9

summary(dx.aov)

p <- plot_richness(ps2.ibs, "BH")
p <- p + geom_boxplot(aes(fill = BH))

dx.aov <- aov(Chao1 ~ BH + SequencingBatch, data1) #0.000111
dx.aov <- aov(Simpson ~ BH + SequencingBatch, data1) #0.32
summary(dx.aov)

summary(lm(Chao1 ~ BH + SequencingBatch, data1) )

# No U
ps2.ibsnoU <- prune_samples(sample_data(ps2.ibs)$BH != "U", ps2.ibs) 
alpha.diversity <- estimate_richness (ps2.ibsnoU, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps2.ibsnoU), alpha.diversity)

p <- plot_richness(ps2.ibsnoU, "BH")
p <- p + geom_boxplot(aes(fill = BH))

dx.aov <- aov(Chao1 ~ BH + SequencingBatch, data1) #0.000111
dx.aov <- aov(Simpson ~ BH + SequencingBatch, data1) #0.32
summary(dx.aov)

summary(lm(Chao1 ~ BH + SequencingBatch, data1) )

# C vd M

ps2.ibsCM <- prune_samples(sample_data(ps2.ibs)$BH == c("C", "M"), ps2.ibs) 
alpha.diversity <- estimate_richness (ps2.ibsCM, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps2.ibsCM), alpha.diversity)

p <- plot_richness(ps2.ibsCM, "BH")
p <- p + geom_boxplot(aes(fill = BH))

dx.aov <- aov(Chao1 ~ BH + SequencingBatch, data1) #0.000111
dx.aov <- aov(Simpson ~ BH + SequencingBatch, data1) #0.32
summary(dx.aov)

summary(lm(Chao1 ~ BH + SequencingBatch, data1) )

ps2.ibsdm <- prune_samples(sample_data(ps2.ibs)$BH == c("C", "M"), ps2.ibs) 
alpha.diversity <- estimate_richness (ps2.ibsdm, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps2.ibsdm), alpha.diversity)

p <- plot_richness(ps2.ibsdm, "BH")
p <- p + geom_boxplot(aes(fill = BH))

dx.aov <- aov(Chao1 ~ BH + SequencingBatch, data1) #0.000111
dx.aov <- aov(Simpson ~ BH + SequencingBatch, data1) #0.32
summary(dx.aov)

# Beta diversity
dist.u <- distance(ps2.ibs, method = "unifrac", weighted = TRUE)
# p <- plot_dist_as_heatmap(dist.bc, title = "unifrac")

ord1 <- ordinate(ps2.ibs, method = "MDS", distance = dist.u)
p <- plot_ordination(ps2.ibs, ord1, color = "Sex")
p <- plot_ordination(ps2.ibs, ord1, color = "SequencingBatch")
p <- plot_ordination(ps2.ibs, ord1, color = "BH")

#dist.b <- distance(ps2.ibs, method = "binary")
clustering <- hclust(dist.u, method = "complete")
#tip.color <- col_factor(palette, levels = levels(Dx))
plot(clustering) #, tip.color = tip.color)


#1. PCoA:  IBS bh
df = as(sample_data(ps2.ibs), "data.frame")
d = distance(ps2.ibs, "unifrac", weighted = FALSE)

B_adonis = adonis(d ~ BH + SequencingBatch, df)
B_adonis # 0.037
p <- plot_ordination(ps2.ibs, ord1, color = "BH")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =4)+stat_ellipse(aes(group=BH))


#0.037

# beta diversity no U

dist.u <- distance(ps2.ibsnoU, method = "unifrac", weighted = TRUE)
# p <- plot_dist_as_heatmap(dist.bc, title = "unifrac")

ord1 <- ordinate(ps2.ibsnoU, method = "MDS", distance = dist.u)
p <- plot_ordination(ps2.ibsnoU , ord1, color = "BH")

#1. PCoA:  IBS and HC
df = as(sample_data(ps2.ibsnoU), "data.frame")
d = distance(ps2.ibsnoU, "unifrac", weighted = FALSE)

B_adonis = adonis(d ~ BH + SequencingBatch, df)
B_adonis # 0.037
p <- plot_ordination(ps2.ibsnoU, ord1, color = "BH")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =4)+stat_ellipse(aes(group=BH))

B_adonis = adonis(d ~ Sex + SequencingBatch, df)
B_adonis
p <- plot_ordination(ps2.ibs, ord1, color = "Sex")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =4)+stat_ellipse(aes(group=Sex))
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Sex              1    0.5232 0.52316  1.8163 0.03905  0.043 *  
# SequencingBatch  2    2.7941 1.39706  4.8502 0.20854  0.001 ***
#  Residuals       35   10.0814 0.28804         0.75242           
#Total           38   13.3987                 1.00000           

#Batches
B_adonis = adonis(d ~ SequencingBatch , df)
B_adonis
p <- plot_ordination(ps2.ibs, ord1, color = "SequencingBatch")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =4)+stat_ellipse(aes(group=SequencingBatch))


# subtypes
hcd = as.dendrogram(clustering)
plot(cut(hcd, h = 0.6)$upper, main = "Upper tree of cut at h=0.6")
groups <- cutree(clustering, k=2)

write.csv(groups, file = "subtypeGroupsIBS.csv")


#3. Association of subtypes with diet 

dietVar <- read.csv("/Users/swapnajoshi/Downloads/lchang_mucosal_microbiome_diet_group_20181127.csv", row.names = 1)

vars <- read.csv("/Users/swapnajoshi/Downloads/sjoshi_microbiome_samples_all_20181127.csv", row.names = 1)

df1 <- merge(sample_data(ps2.ibs), dietVar1, by = "row.names", all.x = TRUE)
row.names(df1) <- df1[,1]
df1$Row.names <- NULL
df1 <- df1[sample_names(ps2.ibs),]
all.equal(row.names(df1), sample_names(ps2.ibs))
sample_data(ps2.ibs) <- df1

# adonis and pcoa with the diet variable
ps3 <- prune_samples(sample_names(ps2.ibs) %in% sample_names(ps2.ibs)[-which(is.na(sample_data(ps2.ibs)$Diet_Group))], ps2.ibs)
sample_data(ps3)$Diet_Group <-  factor(sample_data(ps3)$Diet_Group)
df = as(sample_data(ps3), "data.frame")
d = distance(ps3, "bray")
B_adonis = adonis(d ~ Diet_Group + SequencingBatch, df)
B_adonis #0.48
dist.bc <- distance(ps3, method = "bray")
ord1 <- ordinate(ps3, method = "MDS", distance = dist.bc)
p <- plot_ordination(ps3, ord1, color = "Diet_Group")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =4)+stat_ellipse(aes(group=Diet_Group))

########################################################################
# 1. does the diet differ between subgroups
ps2.ibs <- prune_samples(sample_data(ps2)$Dx == "IBS", ps2) 
diet <- read.csv("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/lchang_mucosal_microbiome_20181126_numericVars.csv", row.names = 1)
diet1 <- diet[row.names(diet)%in%row.names(sample_data(ps2.ibs)),16:30]

df2 <- merge(sample_data(ps2.ibs), diet1, by = "row.names", all.x = TRUE)
row.names(df2) <- df2[,1]
df2$Row.names <- NULL
df2 <- df2[sample_names(ps2.ibs),]
all.equal(row.names(df2), sample_names(ps2.ibs))
sample_data(ps2.ibs) <- df2
if (all.equal(names(groups), sample_names(ps2.ibs))) {sample_data(ps2.ibs)$subtype1 <- groups}

sample_data(ps2.ibs)$subtype1 <- as.factor(sample_data(ps2.ibs)$subtype1)

summary(lm(sample_data(ps2.ibs)$Protein...g ~ sample_data(ps2.ibs)$subtype1, na.action=na.omit))[[4]][2,4]

df4 <- data.frame(sample_data(ps2.ibs))

names1 <- colnames(df4)[6:20]
modelList    <- lapply(names1, function(resp) {
  mF <- formula(paste(resp, " ~ subtype1 + SequencingBatch"))
  summary(lm(mF, data = df4))[[4]][2,4]
})

modelList_df <- as.matrix(modelList)
row.names(modelList_df) <- names1
colnames(modelList_df)  <- "p value"
p value   
Protein...g                                 0.6648971 
Total.fat...g                               0.8771726 
Saturated.fat...g                           0.8827858 
Monounsaturated.fat...g                     0.8519478 
Polyunsaturated.fat...g                     0.3218023 
Carbohydrate...g                            0.4226265 
Total.number.of.grain.servings              0.7305136 
Number.of.whole.grain.servings              0.9959811 
Number.of.non.whole.grain.servings          0.6949362 
Total.number.of.vegetable.servings          0.5536763 
Total.number.of.fruit.servings              0.05269422
Oz.lean.meat.from.beef.pork.lamb.etc        0.9250842 
Oz.lean.meat.from.organ.meats               0.1027332 
Oz.lean.meat.from.franks.luncheon.meats     0.902171  
Total.Oz.lean.meat.from.red.processed.meats 0.9412228   

write.csv(modelList_df, file = "subtype_diet_associations")
df4 $subtype1 <- as.factor(df4$subtype1)
p <- ggplot(aes(x = subtype1, y = Oz.lean.meat.from.organ.meats), data = df4)
p + geom_boxplot(aes(fill = subtype1)) + geom_jitter()

########################################################################
# 2. clinical charachtereistics with subgroups
ps2.ibs <- prune_samples(sample_data(ps2)$Dx == "IBS", ps2) 
clin <- read.csv("/Users/swapnajoshi/Downloads/sjoshi_microbiome_samples_all_20181127.csv", row.names = 1)
clin1 <- clin[row.names(clin)%in%row.names(sample_data(ps2.ibs)),c(17,18,19,21,22,23,25,26,27,39,40,41,42,43,44,45,46,47,48,49,53,55,65,66,70,74,109)]

df2 <- merge(sample_data(ps2.ibs), clin1, by = "row.names", all.x = TRUE)
row.names(df2) <- df2[,1]
df2$Row.names <- NULL
df2 <- df2[sample_names(ps2.ibs),]
all.equal(row.names(df2), sample_names(ps2.ibs))
sample_data(ps2.ibs) <- df2
if (all.equal(names(groups), sample_names(ps2.ibs))) {sample_data(ps2.ibs)$subtype1 <- groups}

sample_data(ps2.ibs)$subtype1 <- as.factor(sample_data(ps2.ibs)$subtype1)

df4 <- data.frame(sample_data(ps2.ibs))

names1 <- colnames(df4)[10:32]
modelList    <- lapply(names1, function(resp) {
  mF <- formula(paste(resp, " ~ subtype1 + SequencingBatch"))
  summary(lm(mF, data = df4))[[4]][2,4]
})

modelList_df <- as.matrix(modelList)
row.names(modelList_df) <- names1
colnames(modelList_df)  <- "p value"
p value  
Age                  0.1035583
BMI                  0.5780575
Education            0.4255301
Marital              0.1289282
Income               0.2228809
BSQ_OverallSx        0.6587186
BSQ_AbdPain          0.6787213
BSQ_Bloating         0.6625906
BSQ_UsualSeverity    0.3554822
BSQ_AgeOnset         0.1600797
BSQdrv_SxDurationYrs 0.5523384
ETI_General_Score    0.379428 
ETI_Physical_Score   0.9715302
ETI_Emotional_Score  0.1720898
ETI_Sexual_Score     0.4210188
ETI_Total_Score      0.8419263
HAD_Anxiety          0.5126736
HAD_Depression       0.688005 
ACE_Score            0.6207461
VSI_Score            0.6393097
PSS_Score            0.7041369
PHQ_Score            0.9669991
IBSSS_Severity       0.412233  

write.csv(modelList_df, file = "subtype_clin_associations")
df4 $subtype1 <- as.factor(df4$subtype1)
p <- ggplot(aes(x = subtype1, y = Age), data = df4)
p + geom_boxplot(aes(fill = subtype1)) + geom_jitter()

p <- ggplot(aes(x = subtype1, y = ETI_Physical_Score), data = df4)
p + geom_boxplot(aes(fill = subtype1)) + geom_jitter()




phylum_glom <- tax_glom(ps2.ibs, taxrank="Phylum")

class_glom <- tax_glom(ps2.ibs, taxrank="Class")

order_glom <- tax_glom(ps2.ibs, taxrank="Order")

family_glom <- tax_glom(ps2.ibs, taxrank="Family")

genus_glom <- tax_glom(ps2.ibs, taxrank="Genus")

species_glom <- tax_glom(ps2.ibs, taxrank="Species")


library("DESeq2")

dds <- phyloseq_to_deseq2(ps2.ibs, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps2.ibs))
write.csv(res1, file = "BH_OTU.csv")

##############
dds <- phyloseq_to_deseq2(phylum_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(phylum_glom), nrow = dim(tax_table(phylum_glom))[1], ncol = dim(tax_table(phylum_glom))[2])))
colnames(res1)[7:13] <- colnames( tax_table(phylum_glom))
colnames(res1)[7] <- "Phylum"
write.csv(res1, "phylumBH.csv")

##############
dds <- phyloseq_to_deseq2(class_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(class_glom), nrow = dim(tax_table(class_glom))[1], ncol = dim(tax_table(class_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(class_glom))
write.csv(res1, "classBH.csv")

##############
dds <- phyloseq_to_deseq2(order_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(order_glom), nrow = dim(tax_table(order_glom))[1], ncol = dim(tax_table(order_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(order_glom))
write.csv(res1, "orderBH.csv")

##################################
dds <- phyloseq_to_deseq2(family_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(family_glom), nrow = dim(tax_table(family_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(family_glom))
write.csv(res1, "family_glomBH.csv")

#######################################
dds <- phyloseq_to_deseq2(genus_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(genus_glom), nrow = dim(tax_table(genus_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(genus_glom))
write.csv(res1, "genus_glomBH.csv")

########################################
dds <- phyloseq_to_deseq2(species_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(species_glom), nrow = dim(tax_table(species_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(species_glom))
write.csv(res1, "species_glomBH.csv")
###############################################################
###############################################################

# CD

ps2.ibscd <- prune_samples(sample_data(ps2.ibs)$BH == c("C","D"), ps2.ibs) 

phylum_glom <- tax_glom(ps2.ibscd, taxrank="Phylum")

class_glom <- tax_glom(ps2.ibscd, taxrank="Class")

order_glom <- tax_glom(ps2.ibscd, taxrank="Order")

family_glom <- tax_glom(ps2.ibscd, taxrank="Family")

genus_glom <- tax_glom(ps2.ibscd, taxrank="Genus")

species_glom <- tax_glom(ps2.ibscd, taxrank="Species")

library("DESeq2")

dds <- phyloseq_to_deseq2(ps2.ibscd, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps2.ibscd))
write.csv(res1, file = "BH_DC_OTU.csv")

##############
dds <- phyloseq_to_deseq2(phylum_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(phylum_glom), nrow = dim(tax_table(phylum_glom))[1], ncol = dim(tax_table(phylum_glom))[2])))
colnames(res1)[7:13] <- colnames( tax_table(phylum_glom))
colnames(res1)[7] <- "Phylum"
write.csv(res1, "phylumBH_CD.csv")

##############
dds <- phyloseq_to_deseq2(class_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(class_glom), nrow = dim(tax_table(class_glom))[1], ncol = dim(tax_table(class_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(class_glom))
write.csv(res1, "classBH_CD.csv")

##############
dds <- phyloseq_to_deseq2(order_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(order_glom), nrow = dim(tax_table(order_glom))[1], ncol = dim(tax_table(order_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(order_glom))
write.csv(res1, "orderBH_CD.csv")

##################################
dds <- phyloseq_to_deseq2(family_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(family_glom), nrow = dim(tax_table(family_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(family_glom))
write.csv(res1, "family_glomBH_CD.csv")

#######################################
dds <- phyloseq_to_deseq2(genus_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(genus_glom), nrow = dim(tax_table(genus_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(genus_glom))
write.csv(res1, "genus_glomBH_CD.csv")

########################################
dds <- phyloseq_to_deseq2(species_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(species_glom), nrow = dim(tax_table(species_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(species_glom))
write.csv(res1, "species_glomBH_CD.csv")

# CM

ps2.ibscm <- prune_samples(sample_data(ps2.ibs)$BH == c("C","M"), ps2.ibs) 

phylum_glom <- tax_glom(ps2.ibscm, taxrank="Phylum")

class_glom <- tax_glom(ps2.ibscm, taxrank="Class")

order_glom <- tax_glom(ps2.ibscm, taxrank="Order")

family_glom <- tax_glom(ps2.ibscm, taxrank="Family")

genus_glom <- tax_glom(ps2.ibscm, taxrank="Genus")

species_glom <- tax_glom(ps2.ibscm, taxrank="Species")

library("DESeq2")

dds <- phyloseq_to_deseq2(ps2.ibscm, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps2.ibscm))
write.csv(res1, file = "BH_DC_OTU.csv")

##############
dds <- phyloseq_to_deseq2(phylum_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(phylum_glom), nrow = dim(tax_table(phylum_glom))[1], ncol = dim(tax_table(phylum_glom))[2])))
colnames(res1)[7:13] <- colnames( tax_table(phylum_glom))
colnames(res1)[7] <- "Phylum"
write.csv(res1, "phylumBH_cm.csv")

##############
dds <- phyloseq_to_deseq2(class_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(class_glom), nrow = dim(tax_table(class_glom))[1], ncol = dim(tax_table(class_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(class_glom))
write.csv(res1, "classBH_cm.csv")

##############
dds <- phyloseq_to_deseq2(order_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(order_glom), nrow = dim(tax_table(order_glom))[1], ncol = dim(tax_table(order_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(order_glom))
write.csv(res1, "orderBH_cm.csv")

##################################
dds <- phyloseq_to_deseq2(family_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(family_glom), nrow = dim(tax_table(family_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(family_glom))
write.csv(res1, "family_glomBH_cm.csv")

#######################################
dds <- phyloseq_to_deseq2(genus_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(genus_glom), nrow = dim(tax_table(genus_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(genus_glom))
write.csv(res1, "genus_glomBH_cm.csv")

########################################
dds <- phyloseq_to_deseq2(species_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(species_glom), nrow = dim(tax_table(species_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(species_glom))
write.csv(res1, "species_glomBH_cm.csv")


# DM

ps2.ibsdm <- prune_samples(sample_data(ps2.ibs)$BH == c("D","M"), ps2.ibs) 

phylum_glom <- tax_glom(ps2.ibsdm, taxrank="Phylum")

class_glom <- tax_glom(ps2.ibsdm, taxrank="Class")

order_glom <- tax_glom(ps2.ibsdm, taxrank="Order")

family_glom <- tax_glom(ps2.ibsdm, taxrank="Family")

genus_glom <- tax_glom(ps2.ibsdm, taxrank="Genus")

species_glom <- tax_glom(ps2.ibsdm, taxrank="Species")

library("DESeq2")

dds <- phyloseq_to_deseq2(ps2.ibsdm, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps2.ibsdm))
write.csv(res1, file = "BH_DC_OTU.csv")

##############
dds <- phyloseq_to_deseq2(phylum_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(phylum_glom), nrow = dim(tax_table(phylum_glom))[1], ncol = dim(tax_table(phylum_glom))[2])))
colnames(res1)[7:13] <- colnames( tax_table(phylum_glom))
colnames(res1)[7] <- "Phylum"
write.csv(res1, "phylumBH_dm.csv")

##############
dds <- phyloseq_to_deseq2(class_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(class_glom), nrow = dim(tax_table(class_glom))[1], ncol = dim(tax_table(class_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(class_glom))
write.csv(res1, "classBH_dm.csv")

##############
dds <- phyloseq_to_deseq2(order_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(order_glom), nrow = dim(tax_table(order_glom))[1], ncol = dim(tax_table(order_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(order_glom))
write.csv(res1, "orderBH_dm.csv")

##################################
dds <- phyloseq_to_deseq2(family_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(family_glom), nrow = dim(tax_table(family_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(family_glom))
write.csv(res1, "family_glomBH_dm.csv")

#######################################
dds <- phyloseq_to_deseq2(genus_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(genus_glom), nrow = dim(tax_table(genus_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(genus_glom))
write.csv(res1, "genus_glomBH_dm.csv")

########################################
dds <- phyloseq_to_deseq2(species_glom, ~ BH)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(species_glom), nrow = dim(tax_table(species_glom))[1], ncol = dim(tax_table(family_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(species_glom))
write.csv(res1, "species_glomBH_dm.csv")