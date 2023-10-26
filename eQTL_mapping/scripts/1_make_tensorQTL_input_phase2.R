# Makeing bed format phenotypes for tensorQTL

#Chr start end ID UNR1 UNR2 UNR3 UNR4 
#chr1 173863 173864 ENSG123 -0.50 0.82 -0.71 0.83
#chr1 685395 685396 ENSG456 -1.13 1.18 -0.03 0.11
#chr1 700304 700305 ENSG789 -1.18 1.32 -0.36 1.26

# Libraries
library(data.table)
library(dplyr)

# Input
args=commandArgs(trailingOnly=TRUE)

omictable_file=args[1]
phe_file=args[2]
anno_file=args[3]
xfam_file=args[4]
output_dir=args[5]
setwd(output_dir)

# Create mapping file to match phenotype to genotype
omictable <- fread(omictable_file, data.table = F)
idmap <- omictable %>%
  select(genotype_individual_id = affymetrix_ID, phenotype_individual_id = RNA_ID) %>%
  filter(!is.na(phenotype_individual_id))

# Read in normalised feature counts
phe <- fread(phe_file, data.table = F)
pheT <- phe %>% select(-V1) %>% t %>% data.frame
names(pheT) <- phe$V1
pheT$feature_id <- rownames(pheT)

anno <- fread(anno_file, data.table = F)

bed <- left_join(pheT, anno%>%select(chromosome, TSS, feature_id)) %>%
  mutate(end=TSS)%>% 
  mutate(start=TSS-1)%>%
  select(-TSS)%>%
  select(Chr = chromosome, start, end, ID = feature_id, everything()) %>%
  arrange(Chr, start) %>%
  filter(!is.na(Chr)) %>%
  rename("#Chr" = Chr)

# Rename IDs to match genotype file
namevec <- base::match(names(bed)[5:ncol(bed)], idmap$phenotype_individual_id)
names(bed)[5:ncol(bed)] <-  as.character(idmap$genotype_individual_id[namevec])

# Sort pheno file IDs
sortedids <- names(bed)[-c(1:4)] %>% sort
bed2 <- bed %>%
  select("#Chr", start, end, ID, sortedids)

fwrite(bed, sep = "\t", file = "INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed")

# Output per-chromosome phenotype files for cis and trans analysis
for(i in as.character(1:22)) {
  bedChr <- bed %>%
    rename(Chr = "#Chr") %>%
    filter(Chr == i)%>%
    rename("#Chr" = Chr)
  fwrite(bedChr, sep = "\t", file = paste0("INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr",i,".bed"))
  rm(bedChr)
}

fwrite(bedChr, sep = "\t", file = paste0("INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr23.bed"))
