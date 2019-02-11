#1. PCoA:  IBS and HC
#2. PCoA: without the IBS-U- is this referring to the IBS only ? I don’t think you need to remove IBS-U unless you are specifically evaluating bowel habit variable
#3. Association of subtypes with Dx, BH, Sex, and other clinical charachteristics (IBS + HC; IBS only)
#4. FDR for dietary correlations--- done
#5. Correlation of microbial subtypes (IBS+HC; IBS only) with diet
#6. Correlation of microbiome with the new diet variable (Mediterranean, standard and modified American, Vegetarian/Vegan, gluten free)

#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")


##################################### correcting batch effects
library(sva)
otuTableMod <- ComBat(as.matrix(t(otu_table(ps))), sample_data(ps)$SequencingBatch)

otuTableMod1 <- matrix(NA, ncol = ncol(otuTableMod), nrow = nrow(otuTableMod))
row.names(otuTableMod1) <- row.names(otuTableMod)
colnames(otuTableMod1) <- colnames(otuTableMod)

for (i in 1:nrow(otuTableMod)) {
  for (j in 1:ncol(otuTableMod)) {
    if (otuTableMod[i,j] <= 0) {
      otuTableMod1[i,j] <- 0
    }
      else if (otuTableMod[i,j] > 0) {
        otuTableMod1[i,j] <- round(otuTableMod[i,j],0)
      }
    }
}

# PS object
psNew <- phyloseq(otu_table(t(otuTableMod1), taxa_are_rows=FALSE), 
               sample_data(ps), tax_table(ps), phy_tree(ps))
psNew2 <- prune_samples(sample_names(psNew) != "A5807", psNew)
psNew.prop <- transform_sample_counts(psNew2, function(otu) otu/sum(otu))

EuclDistAT <- dist(otu_table(psNew.prop), method = "binary")
ord.max <- ordinate(psNew.prop, method = "PCoA", distance = EuclDistAT)
plot_ordination(psNew.prop, ord.max, color="SequencingBatch", title="PCA binary")
hClustering <- hclust(EuclDistAT)
plot(hClustering, hang = -1)

ord.nmds.uni <- ordinate(psNew.prop, method = "PCoA", distance = "unifrac",  weighted = FALSE)
hClustering <- hclust(EuclDistAT)
plot_ordination(psNew.prop, ord.nmds.uni, color="SequencingBatch", title="PCA binary")
plot(hClustering, hang = -1)

# Batch effects correction didnot work well : proceed without correcting
##########################################################################
library(ggplot2)
ps2 <- prune_samples(sample_names(ps) != "A5807", ps)

ps.rare <- rarefy_even_depth(ps2, rngseed = 2210)
ps2.prop <- transform_sample_counts(ps.rare, function(otu) otu/sum(otu))

p <- plot_bar (ps2.prop, fill = "Phylum")
p <- p + facet_wrap(~Dx, scales = "free_x", nrow = 1)

p <- plot_richness(ps2, "Dx")
p <- p + geom_boxplot(aes(fill = Dx))

alpha.diversity <- estimate_richness (ps2, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps2), alpha.diversity)

dx.aov <- aov(Observed ~ Dx + SequencingBatch, data1) #0.473
dx.aov <- aov(Chao1 ~ Dx + SequencingBatch, data1) #0.473
dx.aov <- aov(Shannon ~ Dx + SequencingBatch, data1) #0.73
dx.aov <- aov(Simpson ~ Dx + SequencingBatch, data1) #0.975

p <- plot_richness(ps2, "Sex")
p <- p + geom_boxplot(aes(fill = Sex))

dx.aov <- aov(Chao1 ~ Sex + SequencingBatch, data1) #0.608
dx.aov <- aov(Simpson ~ Sex + SequencingBatch, data1) #0.826

summary(dx.aov)
#Df Sum Sq Mean Sq F value Pr(>F)
#Dx            1    241   240.5   0.194   0.66
#Residuals   154 190936  1239.8      

# bowel habit subtype

p <- plot_richness(ps2, "BH")
p <- p + geom_boxplot(aes(fill = BH))

dx.aov <- aov(Chao1 ~ BH + SequencingBatch, data1) #0.608
dx.aov <- aov(Simpson ~ BH + SequencingBatch, data1) #0.826
summary(dx.aov)


# Beta diversity
dist.u <- distance(ps2, method = "unifrac", weighted = TRUE)
# p <- plot_dist_as_heatmap(dist.bc, title = "unifrac")

ord1 <- ordinate(ps2, method = "MDS", distance = dist.u)
p <- plot_ordination(ps2, ord1, color = "Dx")
p <- plot_ordination(ps2, ord1, color = "SequencingBatch")
p <- plot_ordination(ps2, ord1, color = "Sex")

#dist.b <- distance(ps2, method = "binary")
clustering <- hclust(dist.u, method = "complete")
#tip.color <- col_factor(palette, levels = levels(Dx))
plot(clustering) #, tip.color = tip.color)


#1. PCoA:  IBS and HC
df = as(sample_data(ps2), "data.frame")
d = distance(ps2, "unifrac", weighted = FALSE)

B_adonis = adonis(d ~ Dx + SequencingBatch, df)
B_adonis
p <- plot_ordination(ps2, ord1, color = "Dx")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =3)+stat_ellipse(aes(group=Dx))
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Dx                1     0.378  0.3784  1.3268 0.00697  0.166    
#SequencingBatch   2    10.527  5.2634 18.4543 0.19402  0.001 ***
 # Residuals       152    43.352  0.2852         0.79901           
#Total           155    54.258                 1.00000  

# p = 0.166

B_adonis = adonis(d ~ Sex + SequencingBatch, df)
B_adonis
p <- plot_ordination(ps2, ord1, color = "Sex")
p + theme_bw() + theme(text = element_text(size = 16)) +
geom_point(size =3)+stat_ellipse(aes(group=Sex))
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Sex              1    0.5232 0.52316  1.8163 0.03905  0.043 *  
 # SequencingBatch  2    2.7941 1.39706  4.8502 0.20854  0.001 ***
#  Residuals       35   10.0814 0.28804         0.75242           
#Total           38   13.3987                 1.00000           

#Batches
B_adonis = adonis(d ~ SequencingBatch , df)
B_adonis
p <- plot_ordination(ps2, ord1, color = "SequencingBatch")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =3)+stat_ellipse(aes(group=SequencingBatch))


# subtypes
hcd = as.dendrogram(clustering)
plot(cut(hcd, h = 0.65)$upper, main = "Upper tree of cut at h=0.65")
groups <- cutree(clustering, k=3)

write.csv(groups, file = "subtypeGroups.csv")


#3. Association of subtypes with diet 

dietVar <- read.csv("/Users/swapnajoshi/Downloads/lchang_mucosal_microbiome_diet_group_20181127.csv", row.names = 1)
vars <- read.csv("/Users/swapnajoshi/Downloads/sjoshi_microbiome_samples_all_20181127.csv", row.names = 1)
dietVar  <-  subset(dietVar, dietVar$Diet_Group != 5)
df1 <- merge(sample_data(ps2), dietVar, by = "row.names", all.x = TRUE)
row.names(df1) <- df1[,1]
df1$Row.names <- NULL
df1 <- df1[sample_names(ps2),]
all.equal(row.names(df1), sample_names(ps2))
sample_data(ps2) <- df1

# adonis and pcoa with the diet variable
ps3 <- prune_samples(sample_names(ps2) %in% sample_names(ps2)[-which(is.na(sample_data(ps2)$Diet_Group))], ps2)
sample_data(ps3)$Diet_Group <-  factor(sample_data(ps3)$Diet_Group)
df = as(sample_data(ps3), "data.frame")
d = distance(ps3, "bray")
B_adonis = adonis(d ~ Diet_Group + SequencingBatch, df)
B_adonis #0.36
dist.bc <- distance(ps3, method = "bray")
ord1 <- ordinate(ps3, method = "MDS", distance = dist.bc)
p <- plot_ordination(ps3, ord1, color = "Diet_Group")
p + theme_bw() + theme(text = element_text(size = 16)) +
geom_point(size =3)+stat_ellipse(aes(group=Diet_Group))

########################################################################
# 1. does the diet differ between subgroups
ps2 <- prune_samples(sample_names(ps) != "A5807", ps)
diet <- read.csv("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/lchang_mucosal_microbiome_20181126_numericVars.csv", row.names = 1)
diet1 <- diet[row.names(diet)%in%row.names(sample_data(ps2)),16:30]

df2 <- merge(sample_data(ps2), diet1, by = "row.names", all.x = TRUE)
row.names(df2) <- df2[,1]
df2$Row.names <- NULL
df2 <- df2[sample_names(ps2),]
all.equal(row.names(df2), sample_names(ps2))
sample_data(ps2) <- df2
if (all.equal(names(groups), sample_names(ps2))) {sample_data(ps2)$subtype1 <- groups}

sample_data(ps2)$subtype1 <- as.factor(sample_data(ps2)$subtype1)

summary(lm(sample_data(ps2)$Protein...g ~ sample_data(ps2)$subtype1, na.action=na.omit))[[4]][2,4]

df4 <- data.frame(sample_data(ps2))


names1 <- colnames(df4)[6:20]
modelList    <- lapply(names1, function(resp) {
  mF <- formula(paste(resp, " ~ subtype1 + SequencingBatch"))
  summary(lm(mF, data = df4))[[4]][2,4]
})

modelList_df <- as.matrix(modelList)
row.names(modelList_df) <- names1
colnames(modelList_df)  <- "p value"

write.csv(modelList_df, file = "subtype_diet_associations")
df4 $subtype1 <- as.factor(df4$subtype1)
p <- ggplot(aes(x = subtype1, y = Oz.lean.meat.from.organ.meats), data = df4)
p + geom_boxplot(aes(fill = subtype1)) + geom_jitter()

########################################################################
# 2. clinical charachtereistics with subgroups
ps2 <- prune_samples(sample_names(ps) != "A5807", ps)
clin <- read.csv("/Users/swapnajoshi/Downloads/sjoshi_microbiome_samples_all_20181127.csv", row.names = 1)
clin1 <- clin[row.names(clin)%in%row.names(sample_data(ps2)),c(17,18,19,21,22,23,25,26,27,39,40,41,42,43,44,45,46,47,48,49,53,55,65,66,70,74,109)]

df2 <- merge(sample_data(ps2), clin, by = "row.names", all.x = TRUE)
row.names(df2) <- df2[,1]
df2$Row.names <- NULL
df2 <- df2[sample_names(ps2),]
all.equal(row.names(df2), sample_names(ps2))
sample_data(ps2) <- df2
if (all.equal(names(groups), sample_names(ps2))) {sample_data(ps2)$subtype1 <- groups}

sample_data(ps2)$subtype1 <- as.factor(sample_data(ps2)$subtype1)

df4 <- data.frame(sample_data(ps2))

names1 <- colnames(df4)[10:32]
modelList    <- lapply(names1, function(resp) {
  mF <- formula(paste(resp, " ~ subtype1 + SequencingBatch"))
  summary(lm(mF, data = df4))[[4]][2,4]
})

modelList_df <- as.matrix(modelList)
row.names(modelList_df) <- names1
colnames(modelList_df)  <- "p value"
p value  
Age                  0.4233403 
BMI                  0.2465354 
Education            0.03109307
Marital              0.6026352 
Income               0.529485  
BSQ_OverallSx        0.8544798 
BSQ_AbdPain          0.795998  
BSQ_Bloating         0.4420631 
BSQ_UsualSeverity    0.1895133 
BSQ_AgeOnset         0.5038498 
BSQdrv_SxDurationYrs 0.1080565 
ETI_General_Score    0.9650214 
ETI_Physical_Score   0.04169813
ETI_Emotional_Score  0.7897761 
ETI_Sexual_Score     0.1823637 
ETI_Total_Score      0.3570154 
HAD_Anxiety          0.7348174 
HAD_Depression       0.8835598 
ACE_Score            0.4623743 
VSI_Score            0.6976618 
PSS_Score            0.6026961 
PHQ_Score            0.1262552 
IBSSS_Severity       0.7444993 

write.csv(modelList_df, file = "subtype_clin_associations")
df4 $subtype1 <- as.factor(df4$subtype1)
p <- ggplot(aes(x = subtype1, y = Education), data = df4)
p + geom_boxplot(aes(fill = subtype1)) + geom_jitter()

p <- ggplot(aes(x = subtype1, y = ETI_Physical_Score), data = df4)
p + geom_boxplot(aes(fill = subtype1)) + geom_jitter()

df.res[i,1] <- summary(aov(df4[,1] ~ df4$subtype1, na.action=na.omit))[[1]][1,5]
#df.res[i,2] <- summary(aov(sample_data(ps2)[,i] ~ sample_data(ps2)$subtype1, na.action=na.omit))[[1]][1,5]
#df.res[i,3] <- summary(aov(sample_data(ps2)[,i] ~ sample_data(ps2)$subtype1, na.action=na.omit))[[1]][1,5]




