#!/bin/bash
###########################################
## NOTE: The bed file must be compressed and indexed with the commands below:
############################################

output_dir=$1

module load ceuadmin/tabix/0.2.6
cd $output_dir

bgzip INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed && tabix -p bed INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz

# Compress and index per-chromosome files for conditional analysis
for i in {1..23}
do
  bgzip INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr${i}.bed && tabix -p bed INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr${i}.bed.gz
done