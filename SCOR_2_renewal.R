#1. PCoA:  IBS and HC
#2. PCoA: without the IBS-U- is this referring to the IBS only ? I don’t think you need to remove IBS-U unless you are specifically evaluating bowel habit variable
#3. Association of subtypes with Dx, BH, Sex, and other clinical charachteristics (IBS + HC; IBS only)
#4. FDR for dietary correlations--- done
#5. Correlation of microbial subtypes (IBS+HC; IBS only) with diet
#6. Correlation of microbiome with the new diet variable (Mediterranean, standard and modified American, Vegetarian/Vegan, gluten free)

#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")

##################################### data import
library(phyloseq)
library(dada2)

# source("https://bioconductor.org/biocLite.R")
# biocLite("igraph")
sampleDat1 <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/score_ibs_all_20181213.csv", row.names = 1)
sampleDat2 <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/score_hc_all_20181213.csv", row.names = 1)
sampleDat11 <- sampleDat1[,colnames(sampleDat1)%in%colnames(sampleDat2)] 
sampleDat22 <- sampleDat2[,colnames(sampleDat2)%in%colnames(sampleDat11)] 
all.equal(colnames(sampleDat11), colnames(sampleDat22))
sampleDatAll <- rbind(sampleDat11, sampleDat22)
save(sampleDatAll, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/sampleDatAll.rda")

microb1 <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/naliboff_scor_cluster_w_constip_r1_16s_20181206.csv", row.names = 1)
microb1 <- subset(microb1, microb1$stool_16s == 1); dim(microb1)
sampleDat1 <- sampleDat1[row.names(sampleDat1)%in%row.names(microb1), ]
sampleDat1 <- sampleDat1[row.names(microb1),]
all.equal(row.names(microb1),row.names(sampleDat1))
microb1 <- microb1[,-c(1:16)]
colnames(microb1) <- substr(colnames(microb1),2,100)

microb2 <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/score_hc_all_16s_20181206.csv", row.names = 1)
microb2 <- subset(microb2, microb2$stool_16s == 1); dim(microb2)
sampleDat2 <- sampleDat2[row.names(sampleDat2)%in%row.names(microb2), ]
sampleDat2 <- sampleDat2[row.names(microb2),]
all.equal(row.names(microb2),row.names(sampleDat2))
microb2 <- microb2[,-c(1:51)]
colnames(microb2) <- substr(colnames(microb2),2,100)

com_otus <- colnames(microb1[,colnames(microb1) %in% colnames(microb2)])

microb1 <- microb1[,colnames(microb1) %in% com_otus]
microb2 <- microb2[,colnames(microb2) %in% com_otus]
microb1 <- microb1[,colnames(microb2)]
all.equal(colnames(microb1),colnames(microb2))
all.equal(colnames(sampleDat1),colnames(sampleDat2))
microb <- rbind(microb1,microb2)
sampleDat <- rbind(sampleDat1, sampleDat2)
sampleDat$OCP <- factor(ifelse(sampleDat$HrmSup == "1", 1, 0))
sampleDat$MensP_Stool <- factor(sampleDat$MensP_Stool)

# MensS	Menses Stage: 1=Premenopausal; 2=Perimenopausal; 3=Postmenopausal	
# MensP s
save(sampleDat, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/sampleDat.rda")

taxa <-  read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/16s_OTUID_Tax_Table.csv", row.names = 1); dim(microb)
taxa <- taxa[row.names(taxa)%in%colnames(microb),]
taxa <- taxa[colnames(microb),]
taxa <- as.matrix(taxa)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
all.equal(colnames(microb), row.names(taxa))

tree1 <- read_tree_greengenes("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/IBSmicrobiomeProject/IBSmicrobiomeAnalysisR/data/technicallyCorrectData/fecalMicrobiome/97_otus.tree")

psNew <- phyloseq(otu_table(microb, taxa_are_rows=FALSE), 
                sample_data(sampleDat), tax_table(taxa), phy_tree(tree1))

ps_PrePost <- prune_samples(sample_data(psNew)$MensS==1 | sample_data(psNew)$MensS==3, psNew)
ps_MensPhase <- prune_samples(sample_data(psNew)$MensP_Stool==2 | sample_data(psNew)$MensP_Stool==3, psNew)

sample_data(ps_PrePost)$MensStage <- factor(ifelse(sample_data(ps_PrePost)$MensS == 1, "PreMenopausal", "PostMenopausal"))
sample_data(ps_MensPhase)$MensPhase <- factor(ifelse(sample_data(ps_MensPhase)$MensP_Stool == 2, "Follicular", "Luteal"))
sample_data_MensPhase <- sample_data(ps_MensPhase)
save(psNew, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/psNew_16sSCOR_dat.rda")
save(ps_MensPhase, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/ps_MensPhase.rda")
save(ps_PrePost , file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/ps_prepost.rda")
save(sample_data_MensPhase,  file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/sample_data_MensPhase.rda")
#psNew2 <- prune_samples(sample_names(psNew) != "A5807", psNew)

psNew.prop <- transform_sample_counts(ps_PrePost, function(otu) otu/sum(otu))

EuclDistAT <- dist(otu_table(psNew.prop), method = "binary")
ord.max <- ordinate(psNew.prop, method = "PCoA", distance = EuclDistAT)
plot_ordination(psNew.prop, ord.max, color="OCP", title="PCA binary")
hClustering <- hclust(EuclDistAT)
plot(hClustering, hang = -1)


############################################################

library(ggplot2)

ps.rare <- rarefy_even_depth(ps_PrePost, rngseed = 2210)
ps2.prop <- transform_sample_counts(ps.rare, function(otu) otu/sum(otu))

p <- plot_bar (ps2.prop, fill = "Phylum")
p <- p + facet_wrap(~MensS, scales = "free_x", nrow = 1)

p <- plot_richness(ps_PrePost, "MensS")
p <- p + geom_boxplot(aes(fill = MensS))

sample_data(ps_PrePost)
alpha.diversity <- estimate_richness (ps_PrePost, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps_PrePost), alpha.diversity)
data1$Group <- as.factor(data1$Group)

dx.aov <- aov(Observed ~ MensStage*Group + OCP, data1) #0.108
summary(dx.aov)
dx.aov <- aov(Chao1 ~ MensStage*Group + OCP, data1) #
summary(dx.aov)
dx.aov <- aov(Shannon ~ MensStage*Group + OCP, data1) #0.05
summary(dx.aov)
dx.aov <- aov(Simpson ~ MensStage*Group + OCP, data1) #0.03
summary(dx.aov)

sample_data(ps_PrePost)$Dx <- as.factor(ifelse(sample_data(ps_PrePost)$Group=="1", "HC", "IBS"))

# Beta diversity
dist.u <- distance(ps_PrePost, method = "unifrac", weighted = TRUE)
# p <- plot_dist_as_heatmap(dist.bc, title = "unifrac")

ord1 <- ordinate(ps_PrePost, method = "MDS", distance = dist.u)
p <- plot_ordination(ps_PrePost, ord1, color = "MensStage")
p
p <- plot_ordination(ps_PrePost, ord1, color = "OCP")
p
p <- plot_ordination(ps_PrePost, ord1, color = "Dx")
p

#dist.b <- distance(ps2, method = "binary")
clustering <- hclust(dist.u, method = "complete")
#tip.color <- col_factor(palette, levels = levels(Dx))
plot(clustering) #, tip.color = tip.color)


#1. PCoA:  IBS and HC
# remove NAs

ps_PrePost1 <- prune_samples(sample_names(ps_PrePost) != "A6136", ps_PrePost)
df <- as(sample_data(ps_PrePost1), "data.frame")
d <- distance(ps_PrePost1, "unifrac", weighted = FALSE)

B_adonis = adonis(d ~ MensStage*Dx + OCP, df)
B_adonis
sample_data(ps_PrePost1)$MensStage <- factor(sample_data(ps_PrePost1)$MensStage)
p <- plot_ordination(ps_PrePost1, ord1, color = "MensStage")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =3)+stat_ellipse(aes(group=MensStage))

# Terms added sequentially (first to last)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# MensStage      1    0.1989 0.19893 1.27428 0.01247  0.065 .
# Dx             1    0.1827 0.18273 1.17051 0.01145  0.128  
# OCP            1    0.1422 0.14216 0.91061 0.00891  0.689  
# MensStage:Dx   1    0.1309 0.13092 0.83866 0.00821  0.794  
# Residuals     98   15.2988 0.15611         0.95896         
# Total        102   15.9535                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#########################################################
# 2. clinical charachteristics 

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

#############################################

phylum_glom <- tax_glom(ps_PrePost1, taxrank="Phylum")

class_glom <- tax_glom(ps_PrePost1, taxrank="Class")

order_glom <- tax_glom(ps_PrePost1, taxrank="Order")

family_glom <- tax_glom(ps_PrePost1, taxrank="Family")

genus_glom <- tax_glom(ps_PrePost1, taxrank="Genus")

species_glom <- tax_glom(ps_PrePost1, taxrank="Species")


library("DESeq2")

dds <- phyloseq_to_deseq2(ps_PrePost1, ~ as.factor(MensStage)*as.factor(Group) + as.factor(OCP) + BMI)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_PrePost1))
write.csv(res1, file = "ps_PrePost_OTU.csv")

############################################## 
library(ggplot2)

ps_PrePost_OCPn <- prune_samples(sample_data(ps_PrePost)$OCP == 0, ps_PrePost)

ps.rare <- rarefy_even_depth(ps_PrePost_OCPn, rngseed = 2210)
ps2.prop <- transform_sample_counts(ps.rare, function(otu) otu/sum(otu))

p <- plot_bar (ps2.prop, fill = "Phylum")
p <- p + facet_wrap(~MensStage, scales = "free_x", nrow = 1)

p <- plot_richness(ps_PrePost_OCPn, "MensS")
p <- p + geom_boxplot(aes(fill = MensS))

# sample_data(ps_PrePost_OCPn)
alpha.diversity <- estimate_richness (ps_PrePost_OCPn, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps_PrePost_OCPn), alpha.diversity)
data1$Group <- as.factor(data1$Group)

dx.aov <- aov(Observed ~ MensStage*Group, data1) #0.108
summary(dx.aov)
dx.aov <- aov(Chao1 ~ MensStage*Group, data1) #
summary(dx.aov)
dx.aov <- aov(Shannon ~ MensStage*Group, data1) #0.05
summary(dx.aov)
dx.aov <- aov(Simpson ~ MensStage*Group, data1) #0.03
summary(dx.aov)

sample_data(ps_PrePost_OCPn)$Dx <- as.factor(ifelse(sample_data(ps_PrePost_OCPn)$Group=="1", "HC", "IBS"))

# Beta diversity
dist.u <- distance(ps_PrePost_OCPn, method = "unifrac", weighted = TRUE)
# p <- plot_dist_as_heatmap(dist.bc, title = "unifrac")

ord1 <- ordinate(ps_PrePost_OCPn, method = "MDS", distance = dist.u)
p <- plot_ordination(ps_PrePost_OCPn, ord1, color = "MensStage")
p
# p <- plot_ordination(ps_PrePost_OCPn, ord1, color = "OCP")

p <- plot_ordination(ps_PrePost_OCPn, ord1, color = "Dx")
p

#dist.b <- distance(ps2, method = "binary")
clustering <- hclust(dist.u, method = "complete")
#tip.color <- col_factor(palette, levels = levels(Dx))
plot(clustering) #, tip.color = tip.color)


#1. PCoA:  IBS and HC
# remove NAs

#ps_PrePost_OCPn1 <- prune_samples(sample_names(ps_PrePost_OCPn) != "A6136", ps_PrePost_OCPn)
df <- as(sample_data(ps_PrePost_OCPn), "data.frame")
d <- distance(ps_PrePost_OCPn, "unifrac", weighted = FALSE)

B_adonis = adonis(d ~ MensStage*Dx, df)
B_adonis
sample_data(ps_PrePost_OCPn)$MensStage <- factor(sample_data(ps_PrePost_OCPn)$MensStage)
p <- plot_ordination(ps_PrePost_OCPn, ord1, color = "MensStage")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =3)+stat_ellipse(aes(group=MensStage))

######################

topN = 50
most_abundant_taxa = sort(taxa_sums(ps_PrePost), TRUE)[1:topN]
print(most_abundant_taxa)

set.seed(2210)
ps_PrePostR = rarefy_even_depth(ps_PrePost, sample.size = 500)
ps_PrePost_prop <- filter_taxa(ps_PrePostR, function(x) mean(x) > 1e-1, TRUE)
ps_PrePost_IBS <- prune_samples(sample_data(ps_PrePost_prop)$Dx == "IBS", ps_PrePost_prop)
ps_PrePost_IBS1 <- transform_sample_counts(ps_PrePost_IBS, function(otu) 100* otu/sum(otu))

ps_PrePost_IBS_genus = tax_glom(ps_PrePost_IBS1, "Genus")

plot_bar(ps_PrePost_IBS1, "as.factor(MensStage)", fill = "Genus") + theme(strip.text.x = element_text(size = 8, angle = 90)) + geom_bar(stat="identity")

ps_PrePost_IBS1 <- merge_samples(ps_PrePost_IBS, "MensStage")
ps_PrePost_IBS_prop <- transform_sample_counts(ps_PrePost_IBS1, function(otu) 100* otu/sum(otu))
ps_PrePost_HC  <- prune_samples(sample_data(ps_PrePostR)$Dx == "HC", ps_PrePostR)
ps_PrePost_HC1 <- merge_samples(ps_PrePost_HC, "MensStage")
ps_PrePost_HC_prop <- transform_sample_counts(ps_PrePost_HC1, function(otu) 100* otu/sum(otu))

p <- plot_bar (ps_PrePost_IBS_prop , x = "as.factor(MensStage)",fill = "Phylum")+ geom_bar(stat="identity")
p <- plot_bar (ps_PrePost_HC_prop ,  x = "as.factor(MensStage)",fill = "Phylum")+ geom_bar(stat="identity")

p <- plot_bar (ps_PrePost_IBS_prop , x = "as.factor(MensStage)",fill = "Genus")+ geom_bar(stat="identity")
p <- plot_bar (ps_PrePost_HC_prop ,  x = "as.factor(MensStage)",fill = "Genus")+ geom_bar(stat="identity")


genus_glom <- tax_glom(ps_PrePost, taxrank="Genus")
set.seed(2210)
genus_glomR <- rarefy_even_depth(genus_glom, sample.size = 500)
ps_PrePost_IBS <- prune_samples(sample_data(genus_glomR)$Dx == "IBS", genus_glomR)
ps_PrePost_IBS1 <- merge_samples(ps_PrePost_IBS, "MensS")
ps_PrePost_IBS_prop <- transform_sample_counts(ps_PrePost_IBS1, function(otu) 100* otu/sum(otu))
ps_PrePost_HC <- prune_samples(sample_data(genus_glomR)$Group == "1", genus_glomR)
ps_PrePost_HC1 <- merge_samples(ps_PrePost_HC, "MensS")
ps_PrePost_HC_prop <- transform_sample_counts(ps_PrePost_HC1, function(otu) 100* otu/sum(otu))
p <- plot_bar (ps_PrePost_IBS_prop , x = "as.factor(MensS)",fill = "Genus")+ geom_bar(stat="identity")
p <- plot_bar (ps_PrePost_HC_prop ,  x = "as.factor(MensS)",fill = "Genus")+ geom_bar(stat="identity")
genus_glom <- tax_glom(ps_PrePost, taxrank="Genus")
topN = 50
most_abundant_taxa = names(sort(taxa_sums(genus_glom), TRUE)[1:topN])
print(most_abundant_taxa)
genus_glom1 <-prune_taxa( most_abundant_taxa, genus_glom)
ps_PrePost_IBS <- prune_samples(sample_data(genus_glom1)$Group == 2, genus_glom1)
ps_PrePost_IBS1 <- merge_samples(ps_PrePost_IBS, "MensStage")
ps_PrePost_IBS_prop <- transform_sample_counts(ps_PrePost_IBS1, function(otu) 100* otu/sum(otu))
ps_PrePost_HC <- prune_samples(sample_data(genus_glom1)$Group == 1, genus_glom1) 
ps_PrePost_HC1 <- merge_samples(ps_PrePost_HC, "MensStage")
ps_PrePost_HC_prop <- transform_sample_counts(ps_PrePost_HC1, function(otu) 100* otu/sum(otu))
p <- plot_bar (ps_PrePost_IBS_prop, x = "as.factor(MensS)",fill = "Genus")+ geom_bar(stat="identity")
p <- plot_bar (ps_PrePost_HC_prop, x = "as.factor(MensS)",fill = "Genus")+ geom_bar(stat="identity")

########################### menstrual phase



library(ggplot2)

ps.rare <- rarefy_even_depth(ps_MensPhase, rngseed = 2210)
ps2.prop <- transform_sample_counts(ps.rare, function(otu) otu/sum(otu))

p <- plot_bar (ps2.prop, fill = "Phylum")
p <- p + facet_wrap(~MensP_Stool, scales = "free_x", nrow = 1)
p
p <- plot_richness(ps_MensPhase, "as.factor(MensP_Stool)")
p <- p + geom_boxplot(aes(fill = as.factor(MensP_Stool)))
p
#sample_data(ps_PrePost)
alpha.diversity <- estimate_richness (ps_MensPhase, measures = c("observed", "Chao1", "Shannon", "Simpson"))
data1 <- cbind(sample_data(ps_MensPhase), alpha.diversity)
data1$Group <- as.factor(data1$Group)

dx.aov <- aov(Observed ~ MensP_Stool*Group + OCP, data1) #0.108
summary(dx.aov)
dx.aov <- aov(Chao1 ~ MensP_Stool*Group + OCP, data1) #
summary(dx.aov)
dx.aov <- aov(Shannon ~ MensP_Stool*Group + OCP, data1) #0.05
summary(dx.aov)
dx.aov <- aov(Simpson ~ MensP_Stool*Group + OCP, data1) #0.03
summary(dx.aov)

sample_data(ps_MensPhase)$Dx <- as.factor(ifelse(sample_data(ps_MensPhase)$Group=="1", "HC", "IBS"))

# Beta diversity
dist.u <- distance(ps_MensPhase, method = "unifrac", weighted = TRUE)
# p <- plot_dist_as_heatmap(dist.bc, title = "unifrac")

ord1 <- ordinate(ps_MensPhase, method = "MDS", distance = dist.u)
p <- plot_ordination(ps_MensPhase, ord1, color = "as.factor(MensP_Stool)")
p
p <- plot_ordination(ps_MensPhase, ord1, color = "OCP")
p
p <- plot_ordination(ps_MensPhase, ord1, color = "Dx")
p

#dist.b <- distance(ps2, method = "binary")
clustering <- hclust(dist.u, method = "complete")
#tip.color <- col_factor(palette, levels = levels(Dx))
plot(clustering) #, tip.color = tip.color)

#1. PCoA:  IBS and HC
# remove NAs

ps_MensPhase1 <- prune_samples(sample_names(ps_MensPhase) != "A6136", ps_MensPhase)
df <- as(sample_data(ps_MensPhase), "data.frame")
d <- distance(ps_MensPhase, "unifrac", weighted = FALSE)

B_adonis = adonis(d ~ MensP_Stool*Group + OCP, df)
B_adonis
is.na(sample_data(ps_MensPhase)$MensP_Stool)  <- sample_data(ps_MensPhase)$MensP_Stool
ps_MensPhase$MensP_Stool <- factor(sample_data(ps_MensPhase)$MensP_Stool)
p <- plot_ordination(ps_MensPhase, ord1, color = "MensP_Stool")
p + theme_bw() + theme(text = element_text(size = 16)) +
  geom_point(size =3)+stat_ellipse(aes(group=MensP_Stool))

# Terms added sequentially (first to last)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# MensStage      1    0.1989 0.19893 1.27428 0.01247  0.065 .
# Dx             1    0.1827 0.18273 1.17051 0.01145  0.128  
# OCP            1    0.1422 0.14216 0.91061 0.00891  0.689  
# MensStage:Dx   1    0.1309 0.13092 0.83866 0.00821  0.794  
# Residuals     98   15.2988 0.15611         0.95896         
# Total        102   15.9535                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########################################################################################

# effect size for sample size calculaions

# metabolomics
# sig_ibs_pre <- ibs_preM[,colnames(ibs_preM)%in%row.names(sigMeta)]
# sig_ibs_post <- ibs_postM[,colnames(ibs_postM)%in%row.names(sigMeta)]
# sig_hc_pre <- hc_preM[,colnames(hc_preM)%in%row.names(sigMeta)]

# # all metabolite data IBS vs HCs
# data1 <- rbind(metab.ibs1,metab.hc1)
# sampleDat_metab1 <- subset(sampleDatAll, sampleDatAll$stool_metabolites == 1)
# data1 <- data1[row.names(sampleDat_metab1),29:939]
# all.equal(row.names(data1), row.names(sampleDat_metab1))
# data1$Dx <- as.factor(sampleDat_metab1$Group)
# data1 <- data1[,c(912,1:911)]
# write.csv(data1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_dx.csv")

# pre-post 
data1$Dx <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$MensS == 1 | sampleDat_metab1$MensS == 3)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$MensS <- sampleDat_metab$MensS
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_prepost.csv")

# pre-post IBS
data1$Dx <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$MensS == 1 & sampleDat_metab1$Group == 2 | sampleDat_metab1$MensS == 3 & sampleDat_metab1$Group == 2)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$MensS <- sampleDat_metab$MensS
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_prepost_ibs.csv")

# pre-post HC
data1$Dx <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$MensS == 1 & sampleDat_metab1$Group == 1 | sampleDat_metab1$MensS == 3 & sampleDat_metab1$Group == 1)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$MensS <- sampleDat_metab$MensS
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/all_metabolite_prepost_HCs.csv")


# Follicular vs Luteal IBS
data1$Dx <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$MensP_Stool == 1 & sampleDat_metab1$Group == 2 | sampleDat_metab1$MensP_Stool == 2 & sampleDat_metab1$Group == 2)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$MensP <- sampleDat_metab$MensP_Stool
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_foll_lut_ibs.csv")

# Follicular vs Luteal HCs
data1$Dx <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$MensP_Stool == 1 & sampleDat_metab1$Group == 1 | sampleDat_metab1$MensP_Stool == 2 & sampleDat_metab1$Group == 1)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$MensP_Stool <- sampleDat_metab$MensP_Stool
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_foll_lut_HC.csv")

# Follicular vs Luteal 
data1$Dx <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$MensP_Stool == 1 | sampleDat_metab1$MensP_Stool == 2 )
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$MensP_Stool <- sampleDat_metab$MensP_Stool
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_foll_lut.csv")

# IBS-C vs healthy controls premenopausal
data1$Menopausal_status <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$BH == 1 & sampleDat_metab1$MensS == 1| sampleDat_metab1$BH == 4 & sampleDat_metab1$MensS == 1)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$BH <- sampleDat_metab$BH
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_premeno_IBSC_HC.csv")

# IBS-C vs healthy controls postmenopausal
data1$Menopausal_status <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$BH == 1 & sampleDat_metab1$MensS == 3| sampleDat_metab1$BH == 4 & sampleDat_metab1$MensS == 3)
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$BH <- sampleDat_metab$BH
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_postmeno_IBSC_HC.csv")

# IBS-C vs healthy controls 
data1$Menopausal_status <- NULL
sampleDat_metab <- subset(sampleDat_metab1, sampleDat_metab1$BH == 1 | sampleDat_metab1$BH == 4 )
data2 <- data1[row.names(data1)%in%row.names(sampleDat_metab),]
data2 <- data1[row.names(sampleDat_metab),]; all.equal(row.names(sampleDat_metab), row.names(data2))
data2$BH <- sampleDat_metab$BH
data2 <- data2[,c(912,1:911)]
write.csv(data2, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/metabolite_IBSC_HC.csv")

################ Power for microbiomemicrobiome
library(devtools)
install_github("brendankelly/micropower", host = "https://api.github.com")
library(micropower)

install.packages("HMP" ,repo="http://cran.r-project.org", dep=TRUE)

#### most abundant OTUs 
###################################################
otus <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/OTU_fecat_raw_data.csv", row.names = 1, check.names = FALSE)
otus1 <- otus[,7:16321]
sampledata <- otus[,1:7]
sampledata1 <- sampledata[c(grep("pre", sampledata$Subgroup) , grep("post", sampledata$Subgroup)),]
sampledata1$Subgroup1 <- factor( matrix(unlist(strsplit(as.character(sampledata1$Subgroup), "-")), ncol=2, byrow = TRUE)[,1])

taxa <-  read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/16s_OTUID_Tax_Table.csv", row.names = 1)
taxa1 <- taxa[row.names(taxa)%in%colnames(otus1),]
taxa1 <- taxa1[colnames(otus1),]
all.equal(row.names(taxa1), colnames(otus1))
taxa1 <- as.matrix(taxa1)
colnames(taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ps <- phyloseq(otu_table(otus1, taxa_are_rows=FALSE), 
               sample_data(sampledata1), tax_table(taxa1), phy_tree(tree1))

ps_ibs_prePost <- prune_samples(sample_data(ps)$Group=="2", ps); # 73 samples
##############################################################
ps_PrePost     <- prune_samples(sample_data(psNew)$MensS==1 | sample_data(psNew)$MensS==3, psNew)
ps_IBS_PrePost <- prune_samples(sample_data(psNew)$MensS==1 & sample_data(psNew)$Group=="2"| sample_data(psNew)$MensS==3& sample_data(psNew)$Group=="2", psNew)
ps_HC_PrePost  <- prune_samples(sample_data(psNew)$MensS==1 & sample_data(psNew)$Group=="1"| sample_data(psNew)$MensS==3& sample_data(psNew)$Group=="1", psNew)

ps_MensPhase     <- prune_samples(sample_data(psNew)$MensP_Stool==2 | sample_data(psNew)$MensP_Stool==3, psNew)
ps_IBS_MensPhase <- prune_samples(sample_data(psNew)$MensP_Stool==2 & sample_data(psNew)$Group=="2"| sample_data(psNew)$MensP_Stool==3 & sample_data(psNew)$Group=="2", psNew)
ps_HC_MensPhase  <- prune_samples(sample_data(psNew)$MensP_Stool==2 & sample_data(psNew)$Group=="1"| sample_data(psNew)$MensP_Stool==3 & sample_data(psNew)$Group=="1", psNew)

ps_IBSC_HC       <- prune_samples(sample_data(psNew)$BH==1 | sample_data(psNew)$BH==4, psNew)
ps_pre_IBSC_HC   <- prune_samples(sample_data(psNew)$MensS==1 & sample_data(psNew)$BH=="1"| sample_data(psNew)$MensS==1 & sample_data(psNew)$BH=="4", psNew)
ps_post_IBSC_HC  <- prune_samples(sample_data(psNew)$MensS==3 & sample_data(psNew)$BH=="1"| sample_data(psNew)$MensS==3 & sample_data(psNew)$BH=="4", psNew)


topN = 100
most_abundant_taxa <- names(sort(taxa_sums(ps_PrePost), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_PrePost)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[101,]  <- sample_data(ps_PrePost)$MensS
microb.HMP1        <- microb.HMP1[c(501,1:500),]
row.names(microb.HMP1)[1] <- "Strata"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/PrePostM_hmp.csv")

ps_PrePost.rare <- rarefy_even_depth(ps_PrePost, rngseed = 2210)
x1 <- calcUJstudy(otu_table(ps_PrePost.rare))


power1 <- bootPower(x1, boot_number = 100, subject_group_vector = c(71,32), alpha = 0.05)

most_abundant_taxa <- names(sort(taxa_sums(ps_IBS_PrePost), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_IBS_PrePost)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_IBS_PrePost)$MensS
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "IBSPrePostM"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBSPrePostM_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_HC_PrePost), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_HC_PrePost)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_HC_PrePost)$MensS
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "HCPrePostM"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/HCPrePostM_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_MensPhase), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_MensPhase)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_MensPhase)$MensP_Stool
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "MensPhase"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/MensPhase_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_IBS_MensPhase), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_IBS_MensPhase)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_IBS_MensPhase)$MensP_Stool
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "IBSMensPhase"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBSMensPhase_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_HC_MensPhase), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_HC_MensPhase)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_HC_MensPhase)$MensP_Stool
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "HCMensPhase"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/HCMensPhase_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_IBSC_HC), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_IBSC_HC)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_IBSC_HC)$BH
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "IBSC_HC"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBSC_HC_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_pre_IBSC_HC), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_pre_IBSC_HC)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_pre_IBSC_HC)$BH
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "PreM_IBSC_HC"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/PreM_IBSC_HC_hmp.csv")

most_abundant_taxa <- names(sort(taxa_sums(ps_post_IBSC_HC), TRUE)[1:topN])
microb             <- as.data.frame(t(otu_table(ps_post_IBSC_HC)))
microb.HMP1        <- microb[row.names(microb) %in% most_abundant_taxa,]
microb.HMP1[201,]  <- sample_data(ps_post_IBSC_HC)$BH
microb.HMP1        <- microb.HMP1[c(201,1:200),]
row.names(microb.HMP1)[1] <- "PostM_IBSC_HC"
write.csv(microb.HMP1, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/PostM_IBSC_HC_hmp.csv")

######## microbiome genus level significance figure
library(DESeq2)
genus_glom <- tax_glom(ps_IBS_PrePost, taxrank="Genus")

dds <- phyloseq_to_deseq2(genus_glom, ~ MensS)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(genus_glom), nrow = dim(tax_table(genus_glom))[1], ncol = dim(tax_table(genus_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(genus_glom))
write.csv(res1, "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/genus_PrevsPostIBS.csv")

#ps.rare <- rarefy_even_depth(ps_ibs_prePost, rngseed = 2210)
genus_glom <- tax_glom(ps.rare, taxrank="Genus")

dds <- phyloseq_to_deseq2(genus_glom, ~ Subgroup1)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, as.data.frame(matrix(tax_table(genus_glom), nrow = dim(tax_table(genus_glom))[1], ncol = dim(tax_table(genus_glom))[2])))
colnames(res1)[7:13] <- colnames(tax_table(genus_glom))
write.csv(res1, "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/genus_PrevsPostIBS.csv")

ibs_sig <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBS_significant_prepost.csv")
taxa2 <- as.data.frame(taxa1)
ibs_sig$Phylum <- c("Proteobacteria", "Euryarchaeota", "NA","Firmicutes", "Bacteroidetes", "Firmicutes", "Firmicutes","Actinobacteria", "Proteobacteria", "Proteobacteria", "Firmicutes")

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
ibs_sig$Q_log10 <- ifelse(ibs_sig$log2FC>0, -1*log(ibs_sig$FDR,10), log(ibs_sig$FDR,10))
x = tapply(ibs_sig$Q_log10, ibs_sig$Phylum, function(x) max(x))
x = sort(x, TRUE)
ibs_sig$Phylum = factor(as.character(ibs_sig$Phylum), levels=names(x))
# Genus order
x = tapply(ibs_sig$Q_log10, ibs_sig$Genus, function(x) max(x))
x = sort(x, TRUE)
# ibs_sig$Abundance <- ibs_sig$baseMean
ibs_sig$Genus = factor(as.character(ibs_sig$Genus), levels=names(x))
p <- ggplot(ibs_sig, aes(x=Genus, y=Q_log10, color=Phylum)) + 
  geom_point(size =3.5) +   
  theme(text = element_text(size=20), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))s
ggsave(p, file = "log2Foldchange_phylum_Dx.png")

# Sample size using firmicutes and alpha diversity:
alpha.diversity <- estimate_richness (ps_PrePost, measures =  "Chao1")
data1 <- cbind(sample_data(ps_PrePost), alpha.diversity)
data1$MensS <- as.factor(data1$MensS)
ave_pre  <- mean(subset(data1, data1$MensS == "1")$Chao1)
sd_pre <- sd(subset(data1, data1$MensS == "1")$Chao1)
ave_post <- mean(subset(data1, data1$MensS == "3")$Chao1)
sd_post <- sd(subset(data1, data1$MensS == "3")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post
sps_PrePost

alpha.diversity <- estimate_richness (ps_IBS_PrePost, measures =  "Chao1")
data1 <- cbind(sample_data(ps_IBS_PrePost), alpha.diversity)
data1$MensS <- as.factor(data1$MensS)
ave_pre  <- mean(subset(data1, data1$MensS == "1")$Chao1)
sd_pre <- sd(subset(data1, data1$MensS == "1")$Chao1)
ave_post <- mean(subset(data1, data1$MensS == "3")$Chao1)
sd_post <- sd(subset(data1, data1$MensS == "3")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post

alpha.diversity <- estimate_richness (ps_HC_PrePost, measures =  "Chao1")
data1 <- cbind(sample_data(ps_HC_PrePost), alpha.diversity)
data1$MensS <- as.factor(data1$MensS)
ave_pre  <- mean(subset(data1, data1$MensS == "1")$Chao1)
sd_pre <- sd(subset(data1, data1$MensS == "1")$Chao1)
ave_post <- mean(subset(data1, data1$MensS == "3")$Chao1)
sd_post <- sd(subset(data1, data1$MensS == "3")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post

alpha.diversity <- estimate_richness (ps_MensPhase, measures =  "Chao1")
data1 <- cbind(sample_data(ps_MensPhase), alpha.diversity)
data1$MensP_Stool <- as.factor(data1$MensP_Stool)
ave_pre  <- mean(subset(data1, data1$MensP_Stool == "2")$Chao1)
sd_pre <- sd(subset(data1, data1$MensP_Stool == "2")$Chao1)
ave_post <- mean(subset(data1, data1$MensP_Stool == "3")$Chao1)
sd_post <- sd(subset(data1, data1$MensP_Stool == "3")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post


alpha.diversity <- estimate_richness (ps_IBS_MensPhase, measures =  "Chao1")
data1 <- cbind(sample_data(ps_IBS_MensPhase), alpha.diversity)
data1$MensP_Stool <- as.factor(data1$MensP_Stool)
ave_pre  <- mean(subset(data1, data1$MensP_Stool == "2")$Chao1)
sd_pre <- sd(subset(data1, data1$MensP_Stool == "2")$Chao1)
ave_post <- mean(subset(data1, data1$MensP_Stool == "3")$Chao1)
sd_post <- sd(subset(data1, data1$MensP_Stool == "3")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post

alpha.diversity <- estimate_richness (ps_HC_MensPhase, measures =  "Chao1")
data1 <- cbind(sample_data(ps_HC_MensPhase), alpha.diversity)
data1$MensP_Stool <- as.factor(data1$MensP_Stool)
ave_pre  <- mean(subset(data1, data1$MensP_Stool == "2")$Chao1)
sd_pre <- sd(subset(data1, data1$MensP_Stool == "2")$Chao1)
ave_post <- mean(subset(data1, data1$MensP_Stool == "3")$Chao1)
sd_post <- sd(subset(data1, data1$MensP_Stool == "3")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post

alpha.diversity <- estimate_richness (ps_IBSC_HC, measures =  "Chao1")
data1 <- cbind(sample_data(ps_IBSC_HC), alpha.diversity)
data1$BH <- as.factor(data1$BH)
ave_pre  <- mean(subset(data1, data1$BH == "1")$Chao1)
sd_pre <- sd(subset(data1, data1$BH == "1")$Chao1)
ave_post <- mean(subset(data1, data1$BH == "4")$Chao1)
sd_post <- sd(subset(data1, data1$BH == "4")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post

alpha.diversity <- estimate_richness (ps_pre_IBSC_HC, measures =  "Chao1")
data1 <- cbind(sample_data(ps_pre_IBSC_HC), alpha.diversity)
data1$BH <- as.factor(data1$BH)
ave_pre  <- mean(subset(data1, data1$BH == "1")$Chao1)
sd_pre <- sd(subset(data1, data1$BH == "1")$Chao1)
ave_post <- mean(subset(data1, data1$BH == "4")$Chao1)
sd_post <- sd(subset(data1, data1$BH == "4")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post

alpha.diversity <- estimate_richness (ps_post_IBSC_HC, measures =  "Chao1")
data1 <- cbind(sample_data(ps_post_IBSC_HC), alpha.diversity)
data1$BH <- as.factor(data1$BH)
ave_pre  <- mean(subset(data1, data1$BH == "1")$Chao1)
sd_pre <- sd(subset(data1, data1$BH == "1")$Chao1)
ave_post <- mean(subset(data1, data1$BH == "4")$Chao1)
sd_post <- sd(subset(data1, data1$BH == "4")$Chao1)
(sd_pre+sd_post)/2
ave_pre
ave_post   

# Phylum 

ibs_sig <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBSC_phylum.csv")

# the comparison is reversed so change the FC values

ibs_sig$Q_log10 <- ifelse(ibs_sig$log2FC>0, log(ibs_sig$FDR,10), -1*log(ibs_sig$FDR,10))

x = tapply(ibs_sig$Q_log10, ibs_sig$X, function(x) max(x))
x = sort(x, TRUE)
ibs_sig$Phylum = factor(as.character(ibs_sig$X), levels=names(x))
# Genus order
x = tapply(ibs_sig$Q_log10, ibs_sig$X, function(x) max(x))
x = sort(x, TRUE)
# ibs_sig$Abundance <- ibs_sig$baseMean
ibs_sig$Phylum = factor(as.character(ibs_sig$Phylum), levels=names(x))

p <- ggplot(ibs_sig, aes(x=Phylum, y=Q_log10, color=Phylum)) + 
  geom_point(size = 5.0) +   
  theme(text = element_text(size=20), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(p, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/log2Foldchange_phylum_IBSC1.png")

# genus

ibs_sig <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBSC_genera.csv")
ibs_sig$Genus <- ibs_sig$X
ibs_sig$Q_log10 <- ifelse(ibs_sig$log2FC>0, log(ibs_sig$FDR,10), -1*log(ibs_sig$FDR,10))
x = tapply(ibs_sig$Q_log10, ibs_sig$Phylum, function(x) max(x))
x = sort(x, TRUE)
ibs_sig$Phylum = factor(as.character(ibs_sig$Phylum), levels=names(x))
# Genus order
x = tapply(ibs_sig$Q_log10, ibs_sig$Genus, function(x) max(x))
x = sort(x, TRUE)
# ibs_sig$Abundance <- ibs_sig$baseMean
ibs_sig$Genus = factor(as.character(ibs_sig$Genus), levels=names(x))
p <- ggplot(ibs_sig, aes(x=Genus, y=Q_log10, color=Phylum)) + 
  geom_point(size =4) +   
  theme(text = element_text(size=18), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
ggsave(p, file = "C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/log2Foldchange_genus_IBSC1.png")


###########################################
load("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/ps_prepost.rda")
ibs_sig <- read.csv("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/IBS_significant_prepost.csv", row.names = 1)

genus_glom <- tax_glom(ps_PrePost, taxrank="Genus")
otu <- otu_table(genus_glom)
colnames(otu) <- data.frame(as.matrix(tax_table(genus_glom)))$Genus
genus_prepost <- otu_table(genus_glom)[otu_table(genus_glom) %in% row.names(ibs_sig),]




load("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/ps_MensPhase.rda")

table(sample_data(ps_MensPhase)$MensP_Stool)
fol: 34
lut: 24

table(sample_data(ps_MensPhase)$Group, sample_data(ps_MensPhase)$MensP_Stool)

load("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/ps_prepost.rda")

sigG <- c("Oxalobacter","Eubacterium", "Butyricimonas","Anaerostipes", "Phascolarctobacterium",
"Haemophilus")

genus_glom     <- tax_glom(ps_PrePost, taxrank = "Genus")
df1            <- otu_table(genus_glom)
colnames(df1) <- factor(substr(as.data.frame(tax_table(genus_glom))$Genus,5,100))
df2 <- df1[,colnames(df1)%in%sigG]
metab.ibs <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/naliboff_scor_cluster_w_constip_r1_stool_metabolites_20181205.csv", row.names = 1, check.names = FALSE,stringsAsFactors=FALSE)
sampleDat <- metab.ibs[row.names(metab.ibs)%in%row.names(df2),c(6, 25, 48,51)]
df3 <- df2[row.names(df2)%in%row.names(sampleDat),]
df3 <- df3[row.names(sampleDat),]
all.equal(row.names(sampleDat), row.names(sampleDat))
df4 <- cbind(df3,sampleDat)

cor1 <- rcorr(as.matrix(df4))
adj.p.mat <- apply(cor1$P,2,p.adjust)
corBrainMetab <- cor1$r[1:5,6:8]
pcorBrainMetab <- cor1$P[1:5,6:8]
# adjPcorBrainMetab <- adj.p.mat[1:20,21:22]
# is.na(adjPcorBrainMetab) <- NA
highCor <- ifelse(pcorBrainMetab <= 0.05, 1,0)
highCor<-as.data.frame(highCor)
highCor$sum1<-apply(highCor,1,sum)
highCor1<-highCor[highCor$sum1>0,]
sum2 <- apply(highCor,2,sum, na.rm = TRUE)
sigCor<-sum2[sum2>0]
sigP <- pcorBrainMetab[row.names(pcorBrainMetab)%in%row.names(highCor1), colnames(corBrainMetab)%in%names(sigCor)]
sigR <- corBrainMetab[row.names(corBrainMetab)%in%row.names(sigP), colnames(corBrainMetab)%in%colnames(sigP)]
 # No correlation between symptom severity and 16s genera that were differnt between pre and post ibs women
########################## brain and pre post ibs correlations

ibs_postM <- metab.ibs1[row.names(metab.ibs1) %in% row.names(subset(sampleDatAll, sampleDatAll$MensS=="3")),match.columns.rows(metab.ibs1, sigmetab)]
hc_preM <- metab.hc1[row.names(metab.hc1) %in% row.names(subset(sampleDatAll, sampleDatAll$MensS=="1")),match.columns.rows(metab.ibs1, sigmetab)]
hc_postM <- metab.hc1[row.names(metab.hc1) %in% row.names(subset(sampleDatAll, sampleDatAll$MensS=="3")),match.columns.rows(metab.ibs1, sigmetab)]
ibsPreM16s <- subset(df4, df4$MensS == 1)
ibsPostM16s <- subset(df4, df4$MensS == 3)

brainDat <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/Copy of dataset_for_swapna_SCORE.csv", row.names = 2, check.names = FALSE)
brainVars <- read.csv("C:/Users/swapnajoshi-admin/Dropbox/IBSmicrobiomeAnalysisR/SCOR_2_2018/pre_postMenopausal_IBS/reducedROIS_for_correlations_Swapna_prePostIBS.csv", row.names = 1, check.names = FALSE)
brainVars1 <- brainDat[,colnames(brainDat) %in% row.names(brainVars),]

#brainVars.pre <- brainVars1[match.rows(brainVars1, preM),]
#brainVars.post <- brainVars1[match.rows(brainVars1, postM),]
brain_ibs_pre  <-  brainVars1[match.rows(brainVars1,ibsPreM16s),]; dim(brain_ibs_pre)
brain_ibs_post <-  brainVars1[match.rows(brainVars1,ibsPostM16s),]; dim(brain_ibs_post)
#brain_hc_pre  <-  brainVars1[match.rows(brainVars1,hc_preM),]; dim(brain_hc_pre)
#brain_hc_post <-  brainVars1[match.rows(brainVars1,hc_postM),]; dim(brain_hc_post)

ibs_pre1 <- ibsPreM16s[match.rows(ibsPreM16s,brain_ibs_pre),]; dim(ibs_pre1)
ibs_post1 <- ibsPostM16s[match.rows(ibsPostM16s,brain_ibs_post),]; dim(ibs_post1)
#hc_pre1 <- hc_preM[match.rows(hc_preM,brain_hc_pre),]; dim(hc_pre1)
#hc_post1 <- hc_postM[match.rows(hc_postM,brain_hc_post),]; dim(hc_post1)

all.equal(row.names(brain_ibs_pre), row.names(ibs_pre1))
all.equal(row.names(brain_ibs_post), row.names(ibs_post1))

datPreIBS  <- cbind(ibs_pre1, brain_ibs_pre)
datPostIBS <- cbind(ibs_post1, brain_ibs_post)
#datPreHC  <- cbind(hc_pre1, brain_hc_pre)
#datPostHC <- cbind(hc_post1, brain_hc_post)

# datPostHC <- cbind(hc_post1, brain_hc_post) # no samples in this group

library(Hmisc)

#### pre vs post IBS significant metabolites with significant brain regions
# Pre IBS
cor1 <- rcorr(as.matrix(datPreIBS))
adj.p.mat <- apply(cor1$P,2,p.adjust)
corBrainMetab <- cor1$r[1:5,10:18]
pcorBrainMetab <- cor1$P[1:5,10:18]
# adjPcorBrainMetab <- adj.p.mat[1:20,21:29]
# is.na(adjPcorBrainMetab) <- NA
highCor <- ifelse(pcorBrainMetab <= 0.05, 1,0)
highCor<-as.data.frame(highCor)
highCor$sum1<-apply(highCor,1,sum)
highCor1<-highCor[highCor$sum1>0,]
sum2 <- apply(highCor,2,sum, na.rm = TRUE)
sigCor<-sum2[sum2>0]
sigP <- as.data.frame(pcorBrainMetab[row.names(pcorBrainMetab)%in%row.names(highCor1), colnames(pcorBrainMetab)%in%names(sigCor)])
colnames(sigP) <- "Haemophilus_P.val"
sigP$Haemophilus_R <- corBrainMetab[row.names(corBrainMetab)%in%row.names(highCor1), colnames(corBrainMetab)%in%names(sigCor)]


#                       Haemophilus_P.val Haemophilus_R
# RS__R_InfPrCS_to_R_PosCS        0.032    -0.272
# RS__R_CS_to_R_InfFGOpp          0.028    -0.278
# RS__R_CS_to_L_MPosCgG_S         0.018    -0.300

cor1 <- rcorr(as.matrix(datPostIBS))
adj.p.mat <- apply(cor1$P,2,p.adjust)
corBrainMetab <- cor1$r[1:5,10:18]
pcorBrainMetab <- cor1$P[1:5,10:18]
# adjPcorBrainMetab <- adj.p.mat[1:20,21:29]
# is.na(adjPcorBrainMetab) <- NA
highCor <- ifelse(pcorBrainMetab <= 0.05, 1,0)
highCor<-as.data.frame(highCor)
highCor$sum1<-apply(highCor,1,sum)
highCor1<-highCor[highCor$sum1>0,]
sum2 <- apply(highCor,2,sum, na.rm = TRUE)
sigCor<-sum2[sum2>0]
sigP <- pcorBrainMetab[row.names(pcorBrainMetab)%in%row.names(highCor1), colnames(pcorBrainMetab)%in%names(sigCor)]
sigR <- corBrainMetab[row.names(corBrainMetab)%in%row.names(sigP), colnames(corBrainMetab)%in%colnames(sigP)]

write.csv(sigP, file = "C:/Users/swapnajoshi-admin/Dropbox/SCORE renewal II/PROJECTS/Project 1/Preliminary Data/brain_metabolites_correlation/Pre_post_IBS/PreMnopausalIBS/pcorBrain16sgeneraPreIBS_new.csv")
write.csv(sigR, file = "C:/Users/swapnajoshi-admin/Dropbox/SCORE renewal II/PROJECTS/Project 1/Preliminary Data/brain_metabolites_correlation/Pre_post_IBS/PreMnopausalIBS/corRBrain16sgeneraPreIBS_new.csv")




