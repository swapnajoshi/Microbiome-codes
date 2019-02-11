
source("https://bioconductor.org/biocLite.R")
biocLite("Rcpp")

library(dada2)

path <- "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/IBSmicrobiomeProject/IBSmicrobiomeAnalysisR/data/rawData/microbiome_allBatches/" # CHANGE to the directory containing the fastq files
list.files(path)

# Extract sample names, identify forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Assess quality of samples by position, consider length to truncate at (entered under truncLen); run one line at a time
plotQualityProfile(fnFs)#[1:2])
plotQualityProfile(fnRs)#[1:2])

# Create filtered dataset
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))  # creates output filenames
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))  # creates output filenames
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(148,148),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE)  # maxEE controls the maximum number of expected errors (i.e. read quality), second number is for filtering of reverse reads
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
# seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]  # can use this command to filter out ASVs outside of a desired size range

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)     # remove chimeras
dim(seqtab.nochim)  # shows number of samples x number of non-chimera sequences
sum(seqtab.nochim)/sum(seqtab)  # fraction of sequences that are not chimeras

# Summary of number of sequences across all the steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# Assign taxonomy to ASVs using Silva database
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/tsdon/Desktop/tsdong/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/tsdon/Desktop/tsdong/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


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

#Combined data to phyloseq
library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),phy_tree(fitGTR$tree))
ps


#Exporting tree
library(ape)
tree1 = phy_tree(ps)
ape::write.tree(tree1, "C:/Users/tsdon/Desktop/tsdong/Sample_data/Miseq3_trimmed/tree1.tree")


# Export data so that it can be converted into BIOM
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
output<-cbind(t(seqtab.nochim), taxonomy)
write.table(output, "C:/Users/tsdon/Desktop/tsdong/Sample_data/Miseq3_trimmed/ASV_count_table.txt", sep="\t", col.names=NA)
##### Need to modify .txt file by typing "#OTU" in the upper left box, can then import into QIIME