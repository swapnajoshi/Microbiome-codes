
#source("https://bioconductor.org/biocLite.R")
#biocLite("ape")

library(dada2)
library("DECIPHER")
library("phangorn")
library(phyloseq)
library(ape)

path <- "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_batch1/" # CHANGE to the directory containing the fastq files list.files(path)

# Extract sample names, identify forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Assess quality of samples by position, consider length to truncate at (entered under truncLen); run one line at a time
plotQualityProfile(fnFs)#[1:2]
plotQualityProfile(fnRs)#[1:2]

files1 <- list.files("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_batch1/")
sample.names <- gsub("_R1_001.fastq.gz", "", files1[grep("R1",files1)] )
  
# Create filtered dataset
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))  # creates output filenames
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))  # creates output filenames
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(148,148),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)  
# maxEE controls the maximum number of expected errors (i.e. read quality), second number is for filtering of reverse reads
out

# Generate model for errors in the sequencing run
errF <- learnErrors(filtFs, multithread=TRUE)   # can adjust parameters nbases (default 100 million) to increase amount of sequence data used to assess errors
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)  # visually check if plotted error models represented by black lines fit the observed datapoints (estimates assuming most abundance sequences don't have errors); red line shows estimated error rates from Q scores

memory.limit(100000)

# Identical sequences combined to reduce processing time
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Determine amplicon sequence variants (ASVs) from each read; note: may need to use a different workflow if many files and too little RAM
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]   # number of ASVs on a per sample basis

# Create ASV count table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))  # distribution of lengths of consensus ASVs, should all be near anticipated amplicon length, e.g. 253 for V4
#158  165  168  169  198  207  216  221  250  251  252  253  254 
#1    2    4    8    1    5    1    1    3    5  207 1415   22 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]  # can use this command to filter out ASVs outside of a desired size range

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)     # remove chimeras
dim(seqtab.nochim)  # shows number of samples x number of non-chimera sequences
sum(seqtab.nochim)/sum(seqtab)  # fraction of sequences that are not chimeras
#[1]  0.8361125

# Summary of number of sequences across all the steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write.csv(track, file = "summary_batch1.csv")

# Assign taxonomy to ASVs using Silva database
taxa <- assignTaxonomy(seqtab.nochim, "/Users/swapnajoshi/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/swapnajoshi/Downloads/silva_species_assignment_v132.fa")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa[is.na(taxa)] <- ""
taxonomy <- paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")

output1<-cbind(t(seqtab.nochim), taxonomy)
write.table(output1, "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_batch1/ASV_count_table_batch1.txt", sep="\t", col.names=NA)
##### Need to modify .txt file by typing "#OTU" in the upper left box, can then import into QIIME

######################################## batch 2
path <- "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_resequenced/"

# Extract sample names, identify forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Assess quality of samples by position, consider length to truncate at (entered under truncLen); run one line at a time
plotQualityProfile(fnFs)#[1:2]
plotQualityProfile(fnRs)#[1:2]

files1 <- list.files("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_resequenced/")
sample.names <- gsub("_R1_001.fastq.gz", "", files1[grep("R1",files1)] )

# Create filtered dataset
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))  # creates output filenames
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))  # creates output filenames
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(148,148),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)  
# maxEE controls the maximum number of expected errors (i.e. read quality), second number is for filtering of reverse reads
out
write.csv(out, file= "errors_batch2.csv")

# Generate model for errors in the sequencing run
errF <- learnErrors(filtFs, multithread=TRUE)   # can adjust parameters nbases (default 100 million) to increase amount of sequence data used to assess errors
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)  # visually check if plotted error models represented by black lines fit the observed datapoints (estimates assuming most abundance sequences don't have errors); red line shows estimated error rates from Q scores

memory.limit(100000)

# Identical sequences combined to reduce processing time
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Determine amplicon sequence variants (ASVs) from each read; note: may need to use a different workflow if many files and too little RAM
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]   # number of ASVs on a per sample basis

# Create ASV count table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))  # distribution of lengths of consensus ASVs, should all be near anticipated amplicon length, e.g. 253 for V4
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]  # can use this command to filter out ASVs outside of a desired size range

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)     # remove chimeras
dim(seqtab.nochim)  # shows number of samples x number of non-chimera sequences
sum(seqtab.nochim)/sum(seqtab)  # fraction of sequences that are not chimeras
# [1] 0.6258176

# Summary of number of sequences across all the steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write(track, file = "summary_batch2.csv")

# Assign taxonomy to ASVs using Silva database
taxa <- assignTaxonomy(seqtab.nochim, "/Users/swapnajoshi/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/swapnajoshi/Downloads/silva_species_assignment_v132.fa")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Export data so that it can be converted into BIOM
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
output2<-cbind(t(seqtab.nochim), taxonomy)
write.table(output2, "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_resequenced/ASV_count_table.txt", sep="\t", col.names=NA)
##### Need to modify .txt file by typing "#OTU" in the upper left box, can then import into QIIME

######################################### latest batch

path <- "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_latest_batch/"

# Extract sample names, identify forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Assess quality of samples by position, consider length to truncate at (entered under truncLen); run one line at a time
plotQualityProfile(fnFs)#[1:2]
plotQualityProfile(fnRs)#[1:2]

files1 <- list.files("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_latest_batch/")
sample.names <- gsub("_R1_001.fastq.gz", "", files1[grep("R1",files1)] )

# Create filtered dataset
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))  # creates output filenames
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))  # creates output filenames
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(148,148),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)  
# maxEE controls the maximum number of expected errors (i.e. read quality), second number is for filtering of reverse reads
out
write.csv(out, file = "expected_errors_batch3.csv" )

# Generate model for errors in the sequencing run
errF <- learnErrors(filtFs, multithread=TRUE)   # can adjust parameters nbases (default 100 million) to increase amount of sequence data used to assess errors
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)  # visually check if plotted error models represented by black lines fit the observed datapoints (estimates assuming most abundance sequences don't have errors); red line shows estimated error rates from Q scores

memory.limit(100000)

# Identical sequences combined to reduce processing time
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Determine amplicon sequence variants (ASVs) from each read; note: may need to use a different workflow if many files and too little RAM
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]   # number of ASVs on a per sample basis

# Create ASV count table
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))  # distribution of lengths of consensus ASVs, should all be near anticipated amplicon length, e.g. 253 for V4
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]  # can use this command to filter out ASVs outside of a desired size range

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)     # remove chimeras
dim(seqtab.nochim)  # shows number of samples x number of non-chimera sequences
sum(seqtab.nochim)/sum(seqtab)  # fraction of sequences that are not chimeras

# Summary of number of sequences across all the steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
write.csv(track, file = "summary_batch3.csv")

# Assign taxonomy to ASVs using Silva database
taxa <- assignTaxonomy(seqtab.nochim, "/Users/swapnajoshi/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/swapnajoshi/Downloads/silva_species_assignment_v132.fa")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Export data so that it can be converted into BIOM
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
output3<-cbind(t(seqtab.nochim), taxonomy)
write.table(output, "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_latest_batch/ASV_count_table.txt", sep="\t", col.names=NA)

##### Need to modify .txt file by typing "#OTU" in the upper left box, can then import into QIIME

########################### rename samples based on manifest, retain the best sample in case of duplicates, merge the ASV counts files

bat1 <- read.delim("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_batch1/ASV_count_table_batch1.txt", sep="\t")
man1 <- read.delim("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_latest_batch/Swapna_mapping_batch1.txt")
colnames(bat1) <- substr(colnames(bat1),2,5)
colnames(bat1) <- gsub("_", "", colnames(bat1))
row.names(bat1) <- bat1[,1]
bat1 <- bat1[,-1]
colnames(bat1)[70] <- "taxonomy"
row.names(man1) <- man1$Sample
man2 <- man1[row.names(man1) %in% colnames(bat1),];  dim(man2)
man2 <- man2[colnames(bat1)[1:69],]
all.equal(row.names(man2), colnames(bat1)[1:69])
colnames(bat1)[1:69] <- as.character(man2$X.SampleID)

##
reseq <- read.delim("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_resequenced/ASV_count_table.txt", sep="\t")
manRes <- read.delim("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData//microbiome_latest_batch/Swapna_resequenced_mapping.txt", sep = "\t")
manRes$Sample <- paste(manRes$Plate, manRes$Well, sep = "")
colnames(reseq) <- substr(colnames(reseq),2,5)
colnames(reseq) <- gsub("_", "", colnames(reseq))
row.names(reseq) <- reseq[,1]
reseq1 <- reseq[,-1]
colnames(reseq1)[35] <- "taxonomy"
row.names(manRes) <- manRes$Sample
manRes2 <- manRes[row.names(manRes) %in% colnames(reseq1),];  dim(manRes2)
manRes2 <- manRes2[colnames(reseq1)[1:34],]
all.equal(row.names(manRes2), colnames(reseq1)[1:34])
colnames(reseq1)[1:34] <- substr(as.character(manRes2$X.SampleID),1,5)

bat12uniq <- bat1[,!colnames(bat1)%in%colnames(reseq1)]
bat12uniq$taxonomy <- bat1$taxonomy
row.names(bat12uniq) <- row.names(bat1)
batch12 <- merge(bat12uniq, reseq1, by = "row.names")
row.names(batch12) <- batch12[,1]
batch12$Row.names <- NULL
batch12$taxonomy.x <- NULL
colnames(batch12)[78] <-  "taxonomy"

##
latBat <- read.delim("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/microbiome_latest_batch/ASV_count_table.txt", sep="\t")
colnames(latBat) <- substr(colnames(latBat),1,5)
row.names(latBat) <- latBat[,1]
latBat <- latBat[,-1]
colnames(latBat)[81] <- "taxonomy"

IBS_microb_ASV_counts <- merge(batch12, latBat, by = "row.names")
row.names(IBS_microb_ASV_counts) <- IBS_microb_ASV_counts[,1]
IBS_microb_ASV_counts <- IBS_microb_ASV_counts[,-c(1,79)]
colnames(IBS_microb_ASV_counts)[158] <- "taxonomy"
write.table(IBS_microb_ASV_counts, file = "/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/IBS_microb_ASV_counts.txt", sep = "\t")

# Tree
seqtab.nochim <- t(IBS_microb_ASV_counts[,1:157])
#creating a phylogenetic tree: process takes a long time 
library("DECIPHER")
library("phangorn")
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))
# pml.control(trace = 1) may be useful to check if optim.pml is still doing something or is crashed
#The most time consuming part will be the stochastic rearrangements
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic", ratchet.par = list(iter = 5L, maxit = 5L, prop = 1/3))
#This is only 5 times instead of (a maximum) of 100 rearrangements. Now you could extrapolate how long phangorn would run and try
fitGTR <- optim.pml(fitGTR, rearrangement = "stochastic")
# this could be parallelized at the cost of memory consumption and in the end optimize all the other parameters again (they should not change much)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "NNI", control = pml.control(trace = 1))
detach("package:phangorn", unload=TRUE)

taxa <- assignTaxonomy(seqtab.nochim, "/Users/swapnajoshi/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/swapnajoshi/Downloads/silva_species_assignment_v132.fa")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

metaDat1 <- read.delim("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/IBS_metadata_all_small.txt", row.names = 1)
metaDat1 <- metaDat1[row.names(seqtab.nochim),]

subject <- row.names(seqtab.nochim)
Dx <- metaDat1$Dx
Sex <- metaDat1$Gender
BH <- metaDat1$BH.at.Colon.Exam
SequencingBatch <- as.factor(metaDat1$batch) 
samdf <- data.frame(Subject=subject, Dx=Dx, Sex= Sex, BH= BH, SequencingBatch=SequencingBatch)
rownames(samdf) <- subject

#Combined data to phyloseq

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), tax_table(taxa), phy_tree(fitGTR$tree))

taxa_names(ps) <- paste("OTU-", seq(1:dim(otu_table(ps))[2]), sep="")
ps

plot_richness(ps, x="Dx", measures=c("Shannon", "Simpson"), color="Dx")

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.uni <- ordinate(ps.prop, method = "PCoA", distance = "unifrac", weighted = FALSE)
plot_ordination(ps.prop, ord.nmds.uni, color="Dx", title="PCA Unifrac")
plot_ordination(ps.prop, ord.nmds.uni, color="BH", title="PCA Unifrac")
plot_ordination(ps.prop, ord.nmds.uni, color="Sex", title="PCA Unifrac")
plot_ordination(ps.prop, ord.nmds.uni, color="Sex", title="PCA Unifrac")
plot_ordination(ps.prop, ord.nmds.uni, color="SequencingBatch", title="PCA Unifrac")

phyloseq_adonis = adonis(d ~ Dx , df)
phyloseq_adonis
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#Dx          1    0.1142 0.114158  1.5777 0.01008  0.157
#Residuals 155   11.2155 0.072358         0.98992       
#Total     156   11.3297                  1.00000   

phyloseq_adonis = adonis(d ~ Sex , df)
phyloseq_adonis
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
#Sex         1    0.0784 0.078383  1.0798 0.00692  0.331
#Residuals 155   11.2513 0.072589         0.99308       
#Total     156   11.3297                  1.00000   

phyloseq_adonis = adonis(d ~ BH , df)
phyloseq_adonis

plot_richness(ps, x="BH", measures=c("Shannon", "Simpson"), color="Dx")
plot_net(ps, maxdist=0.4, color="Dx", shape="Dx")
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# BH          4    0.4345 0.108635  1.5156 0.03835  0.047 *
#  Residuals 152   10.8951 0.071679         0.96165         
# Total     156   11.3297                  1.00000    

# Taxa level differences

phylum_glom <- tax_glom(ps, taxrank="Phylum")

class_glom <- tax_glom(ps, taxrank="Class")

order_glom <- tax_glom(ps, taxrank="Order")

family_glom <- tax_glom(ps, taxrank="Family")

genus_glom <- tax_glom(ps, taxrank="Genus")

species_glom <- tax_glom(ps, taxrank="Species")

### Autonomic data 
# pain bloating severity 
# correlate with microbiome.

# Diet and microbiome

library("DESeq2")

dds <- phyloseq_to_deseq2(ps, ~ Dx)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps))
write.csv(res1, file = "Dx_OTU.csv")

##############
dds <- phyloseq_to_deseq2(phylum_glom, ~ Dx)
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
write.csv(res1, "phylumDx.csv")

##############
dds <- phyloseq_to_deseq2(class_glom, ~ Dx)
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
write.csv(res1, "classDx.csv")

##############
dds <- phyloseq_to_deseq2(order_glom, ~ Dx)
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
write.csv(res1, "orderDx.csv")

##################################
dds <- phyloseq_to_deseq2(family_glom, ~ Dx)
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
write.csv(res1, "family_glomDx.csv")

#######################################
dds <- phyloseq_to_deseq2(genus_glom, ~ Dx)
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
write.csv(res1, "genus_glomDx.csv")

########################################
dds <- phyloseq_to_deseq2(species_glom, ~ Dx)
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
write.csv(res1, "species_glomDx.csv")

############################## subtypes
ps.prop2 <- prune_samples(sample_names(ps.prop) != "A5807", ps.prop)

EuclDistAT <- dist(otu_table(ps.prop2), method = "binary")
ord.nmds.uni <- ordinate(ps.prop2, method = "PCoA", distance = "unifrac",  weighted = FALSE)
hClustering <- hclust(EuclDistAT)
plot_ordination(ps.prop, ord.nmds.uni, color="SequencingBatch", title="PCA binary")
plot(hClustering, hang = -1)

hcd = as.dendrogram(hClustering)
plot(cut(hcd, h = 0.95)$upper, main = "Upper tree of cut at h=0.95")
groups <- cutree(hClustering, k=3)
ps

write.csv(groups, file = "subtypeGroups.csv")

############################ subgroup correlations
ps2 <- prune_samples(sample_names(ps) != "A5807", ps)
For Swapna ANS Data.csv

############################ ANS correlons
ans <- read.csv("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/For Swapna ANS Data.csv", row.names = 2)
ans1 <- ans[row.names(ans)%in%row.names(sample_data(ps2)),]
psMod <- prune_samples(sample_names(ps) %in% row.names(ans1), ps)
ans1 <- ans1[sample_names(psMod),]
ps_ans <- phyloseq(otu_table(psMod, taxa_are_rows=FALSE), 
                   sample_data(ans1), tax_table(psMod))

################################ diet correlations
diet <- read.csv("/Users/swapnajoshi/Dropbox/IBSmicrobiomeAnalysisR/data/rawData/lchang_mucosal_microbiome_20181126_numericVars.csv", row.names = 1)
diet1 <- diet[row.names(diet)%in%row.names(sample_data(ps2)),]
psMod <- prune_samples(sample_names(ps) %in% row.names(diet1), ps)
diet1 <- diet1[sample_names(psMod),]
ps_diet <- phyloseq(otu_table(psMod, taxa_are_rows=FALSE), 
                   sample_data(diet1), tax_table(psMod))


r <- cor(otu_table(ps_diet), sample_data(ps_diet), use = "pairwise.complete.obs") # MUCH MUCH faster than corr.test()
cor2pvalue = function(r, n) {
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  se <- sqrt((1-r*r)/(n-2))
  out <- list(r, n, t, p, se)
  names(out) <- c("r", "n", "t", "p", "se")
  return(out)
}

# get a list with matrices of correlation, pvalues, standard error, etc.
result = as.data.frame(cor2pvalue(r,131))
row.names(result)==row.names(tax_table(psMod)) 
result1 <- cbind (result, tax_table(psMod))
adjp <- as.data.frame(matrix(NA, nrow = nrow(result1), ncol = 30))

for (i in 1:30) {
  adjp[,i] <- p.adjust(result1[,61+i], method = "BH")
}

colnames(adjp) <- paste("adj_",colnames(result1)[62:91], sep = "")
row.names(adjp) <- row.names(result1)
result2 <- cbind(result1, adjp)                   
write.csv(result2, file = "dietaryCorrelations_fdr.csv")

##### t-.test with discrete diet variables  [1] "Diet_StandardAmerican"   "Diet_ModAmerican"       
[3] "Diet_Mediterranean"      "Diet_Paleo"             
[5] "Diet_Vegan"              "Diet_Vegetarian"        
[7] "Diet_LactoVegetarian"    "Diet_OvoVegetarian"     
[9] "Diet_LactoOvoVegetarian" "Diet_Pescatarian"       
[11] "Diet_RawVegan"           "Diet_GlutenFree"        
[13] "Diet_DairyFree"          "Diet_FODMAP"   
#1.
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,1])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_StandardAmerican <- as.factor(sample_data(ps_diet_noNA)$Diet_StandardAmerican)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_StandardAmerican)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_StandardAmerican_OTU.csv")

#2. 
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,2])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_ModAmerican <- as.factor(sample_data(ps_diet_noNA)$Diet_ModAmerican)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_ModAmerican)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_ModAmerican_OTU.csv")

# 3.
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,3])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_Mediterranean <- as.factor(sample_data(ps_diet_noNA)$Diet_Mediterranean)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_Mediterranean)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_Mediterranean_OTU.csv")

#4.
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,4])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_Paleo <- as.factor(sample_data(ps_diet_noNA)$Diet_Paleo)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_Paleo)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_Paleo_OTU.csv")

#5.
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,5])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_Vegan <- as.factor(sample_data(ps_diet_noNA)$Diet_Vegan)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_Vegan)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_Vegan_OTU.csv")

#6.
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,6])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_Vegetarian <- as.factor(sample_data(ps_diet_noNA)$Diet_Vegetarian)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_Vegetarian)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_Vegetarian_OTU.csv")

#12
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,12])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_GlutenFree <- as.factor(sample_data(ps_diet_noNA)$Diet_GlutenFree)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_GlutenFree)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_GlutenFree_OTU.csv")

#14
ps_diet_noNA <- prune_samples(sample_names(sample_data(ps_diet)[which(!is.na(sample_data(ps_diet)[,14])),]), ps_diet)
sample_data(ps_diet_noNA)$Diet_FODMAP <- as.factor(sample_data(ps_diet_noNA)$Diet_FODMAP)
dds <- phyloseq_to_deseq2(ps_diet_noNA, ~ Diet_FODMAP)
table(colSums(counts(dds) == 0))
idx <- which.max(colSums(counts(dds) == 0))
dds2 <- dds[ , -idx]

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds2), 1, gm_mean)
diagdds = estimateSizeFactors(dds2, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res1 <- results(diagdds)
res1 <- cbind(res1, tax_table(ps_diet_noNA))
write.csv(res1, file = "Diet_FODMAP_OTU.csv")

################### heatmap
#plot_heatmap(ps, "NMDS", "bray", "Dx", "Family")

## correlation plots
franksDf <- result2[which(result2$adj_p.Oz.lean.meat.from.franks.luncheon.meats<0.05),c(29,157,123,127)]
wholeGrainDf <- result2[which(result2$adj_p.Number.of.whole.grain.servings<0.05),c(23,151,123,127)]
colnames(franksDf) <- c("Correlation_coefficient", "FDR P-value", "Phylum", "Genus")
colnames(wholeGrainDf) <- c("Correlation_coefficient", "FDR P-value", "Phylum", "Genus")
#sigtab = franksDf
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(franksDf, !is.na(Genus))

# Phylum order
x = tapply(sigtabgen$Corrrelation_coefficient, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$Corrrelation_coefficient, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = Correlation_coefficient, color = Phylum)) + geom_point(size=4) + 
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))


sigtabgen = subset(wholeGrainDf, !is.na(Genus))

# Phylum order
x = tapply(sigtabgen$Corrrelation_coefficient, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$Corrrelation_coefficient, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = Corrrelation_coefficient, color = Phylum)) + geom_point(size=4) + 
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))





scale_color_manual(values=c("red", "blue", "green", "magenta", "black", "gray", "orange"))
         


scale_color_manual(values=c("red",  "black", "orange"))

sigtabgen = subset(wholeGrainDf, !is.na(Genus))
                    
# Phylum order
x = tapply(sigtabgen$Correlation_coefficient, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$Correlation_coefficient, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
p <- ggplot(sigtabgen, aes(x = Genus, y = Correlation_coefficient, color = Phylum)) + geom_point(size=4)
p + scale_color_manual(values=c("red", "blue", "green", "magenta", "black", "gray", "orange")) + 
theme(text = element_text(size=16), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
                    
                    
sigtabgen = subset(franksDf, !is.na(Genus))                   
sigtabgen$Phylum <- factor(sigtabgen$Phylum, levels = c("Firmicutes", "Bacteroidetes", "Proteobacteria"))
# Phylum order
x = tapply(sigtabgen$Correlation_coefficient, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$Correlation_coefficient, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
p <- ggplot(sigtabgen, aes(x = Genus, y = Correlation_coefficient, color = Phylum)) + geom_point(size=4)
p + scale_color_manual(values=c("green", "red", "blue")) + 
  theme(text = element_text(size=16), axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))






