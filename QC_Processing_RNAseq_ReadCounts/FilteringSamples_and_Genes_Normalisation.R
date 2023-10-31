#-----------------------------------------------------
# Filtering poor samples, lowly expressed genes 
# and creating normalised data
# TMM normalised FPKM values are further inverse rank normalised
#-----------------------------------------------------

# Load libraries
library(ggplot2)
library(data.table)
library(dplyr)  

# set working directory
setwd("~/06_AllSamples_Final_Analysis/Data/")
#---------------------------------------------------------------------------
# 1. Read the appropriate files
#---------------------------------------------------------------------------
depth_info <- read.csv("Covariates/Processed/raw_and_globin_depth_per_sample.csv")
count <- fread("Counts/Processed/5591-star-fc-genecounts_added_geneinfo.csv", data.table=FALSE)

# get sequencing covariates 
meta <- fread("Covariates/Processed/INTERVAL_RNA_batch1-15_master_covariates_finalrelease_2021_02_19.csv", data.table=FALSE)
meta1 <- merge(meta, depth_info, by.x = "sample_id", by.y ="Sample_id")
rownames(meta1) <- meta1$RNA_id
meta2 <- meta1[colnames(count[12:ncol(count)]),]
table(meta2$RNA_id == colnames(count[12:ncol(count)]))


#---------------------------------------------------------------------------
# 2. For batch 1 and Batch 15 remove the single run (run1/run2) riles 
# keep the combined run
#--------------------------------------------------------------------------
# keep samples with 2x reads and remove 1x samples from batch 1 and batch 15
run1 <- colnames(count)[grepl("_1", colnames(count))]
run2 <- colnames(count)[grepl("_2", colnames(count))]
run_12 <- c(run1, run2)

count_sub <- count[, -which(colnames(count) %in% run_12)]
metadata <- meta2[-which(meta2$RNA_id %in% run_12),]
metadata$RIN <- as.numeric(metadata$Agilent_RINe)

table(metadata$RNA_id == colnames(count_sub[12:ncol(count_sub)]))


#-----------------------------------------------------
# 3. Remove individuals
#-----------------------------------------------------
# There are inidividuals that have NA values for sequencing covariates
# These are the same individuals removed previously owing to batch two
#N=8 samples
NAs_sequenicng_info <- metadata$sample_id[which(is.na(metadata$sequencingBatch))]

#Remove samples with RIN < 4
# There are N=15 individuals
remove_RIN_Less4 <- as.character(metadata$sample_id[which(metadata$RIN < 4)])

#Remove samples with ReadDepth Less than 10
# There are N=10 individuals
remove_ReadDepth_lessthan_10Million <- metadata$sample_id[which(metadata$depth_rawdata < 10000000)]

#Remove N=6 related individuals based on genetic PCs
metadata$affy <- as.character(metadata$affymetrix_ID)
related_affy_data <- read.table("related_samples.txt")
related_affy <- related_affy_data[,1]
remove_related_RNAid <- metadata$sample_id[which(metadata$affy %in% related_affy)]

# Total remove N=39 individuals
remove_ind <- c(NAs_sequenicng_info, remove_RIN_Less4, remove_ReadDepth_lessthan_10Million, remove_related_RNAid)

# Remove individuals 
count_filtSamp <- count_sub[, -which(colnames(count_sub) %in% remove_ind)]

#-----------------------------------------------------
# 3. Remove globin genes % across samples
#-----------------------------------------------------
# remove globin genes
globin_genes <- c("HBA1","HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "MB")
read_depth <- data.frame(colSums(count_filtSamp[,12:ncol(count_filtSamp)]))
colnames(read_depth) <- "depth"
a <- count_filtSamp[which(count_filtSamp$gene_symbol %in% globin_genes),]
rownames(a) <- a$gene_symbol
a1 <- a[, 12:ncol(a)]
counts_globin <- t(a1)
counts_per <- NULL
for(i in 1:ncol(counts_globin)){
  print(i)
  b <- colnames(counts_globin)[i]
  b1 <- paste0(b, "_Per")
  per <- data.frame(counts_globin[,i]/read_depth$depth * 100)
  colnames(per) <- b1
  counts_per <- cbind(counts_per, per[,1])
  colnames(counts_per)[i] <- b1
}
rownames(counts_per) <- rownames(counts_globin)
counts_per <- data.frame(counts_per)

#-----------------------------------------------------
# 4. Remove globin genes, rRNA, and lowly expressed genes
# Create DGEList object y and remove outlier genes
#-----------------------------------------------------
rRNA <- count_filtSamp$gene_symbol[which(grepl("rRNA", count_filtSamp$class))]
rRNA_globin <- c(rRNA, globin_genes)
count_filtSampGenes  <-count_filtSamp[-which(count_filtSamp$gene_symbol %in% rRNA_globin),]

group <- c(rep("A",length(12:ncol(count_filtSampGenes))))

# creat DGEList
y <- edgeR::DGEList(counts=count_filtSampGenes[,12:ncol(count_filtSampGenes)], genes=count_filtSampGenes[,1:11], group = group)

# Filter lowly expressed genes
# Select genes with > 0.5 CPM in at least 10% of the samples
keep <- rowSums(cpm(y) > 0.5) >= round(ncol(y$counts)*0.01) #10% is N=276 individuals

summary(keep)
y1 <- y[keep, , keep.lib.sizes = FALSE]


# TMM Normalisation
y2 <- calcNormFactors(y1, method = "TMM")

#TMM_log2CPM <- cpm(y2, log = TRUE, normalized.lib.sizes = TRUE)
TMM_log2FPKM <- rpkm(y2, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25, gene.length = y2$genes$gene_length)
TMM_log2FPKM1 <- cbind(y2$genes, TMM_log2FPKM)
write.csv(TMM_log2FPKM1, "INTERVAL_FilteredSamplesGenes_TMMNormalised_FPKM_values", row.names = FALSE, quote = FALSE)

#-----------------------------------------------------------
# 5. Swap Column Names
#-----------------------------------------------------------
TMM_log2FPKM1 <- fread("Counts/Processed/INTERVAL_FilteredSamplesGenes_TMMNormalised_FPKM_values.csv", data.table=FALSE)
swap <- read.csv("SamplesSwaps_to_Change_FromMBV.csv")
remove_mismatch <- as.character(swap$BamFile_Old[which(swap$Status == "remove")])
swap1 <- swap[which(swap$Status == "swap"),]
swap1$TempName <- paste0(swap1$BamFile_New, "_a")
  
# remove samples not showing high concordance with their corresponding genotype or any other genotype 
TMM_log2FPKM1_remove <- TMM_log2FPKM1[,-which(colnames(TMM_log2FPKM1) %in% remove_mismatch)]
dat <- TMM_log2FPKM1_remove 
for(i in 1:nrow(swap1)){
  old <- swap1$BamFile_Old[i]
  colnames(dat)[which(colnames(dat) == old)] <- swap1$TempName[i]
}
colnames(dat) <- gsub("_a", "", colnames(dat))

TMM_log2FPKM1_swapsSwapped_mismatchRemoved <- dat
path <- "~/06_AllSamples_Final_Analysis/Data/Counts/Processed/"
write.csv(TMM_log2FPKM1_swapsSwapped_mismatchRemoved, paste0(path,"INTERVAL_FilteredSamplesGenes_TMMNormalisedFPKM_swapsSwapped_mismatchRemoved.csv"), row.names = FALSE, quote = FALSE)


#-----------------------------------------------------------
# 6. Inverse Rank Normalise data 
#-----------------------------------------------------------

## file preprocessing - gene expression files
library(GenABEL)
gene <- TMM_log2FPKM1_swapsSwapped_mismatchRemoved[,12:ncol(TMM_log2FPKM1_swapsSwapped_mismatchRemoved)]


#z-transform the expression values for each gene were transformed into a standard normal based on rank (to minimize the effects of outliers on the regression scores). 
rownames(gene) <- TMM_log2FPKM1_swapsSwapped_mismatchRemoved$gene_id
gene_transpose <- t(gene)
gene_RankTransform <- apply(gene_transpose, 2, rntransform)
rownames(gene_RankTransform) <- rownames(gene_transpose)
write.csv(gene_RankTransform, "GeneExpr_PEER_TmmInvRankNormalised_swapsSwapped_mismatchRemoved.csv", row.names =TRUE, quote = FALSE) # this data was used for PEER factor analysis

#---------------------------------------------------


