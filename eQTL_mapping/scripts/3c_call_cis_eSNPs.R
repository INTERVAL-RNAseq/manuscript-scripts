# Filtering significant cis-egenes and cis-esnps

args=commandArgs(trailingOnly=TRUE)
pergene=args[1]
prefix_nominal=args[2]
prefix_bim=args[3]
prefix_results=args[4]
b37_b38_map_file=args[5]
info_file=args[6]

# pergene="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cisPerGene_ALLchr.csv"
# prefix_nominal="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cisNominal_chr"
# prefix_bim="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/genotypes/INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_chr"
# prefix_results="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cis_chr"
b37_b38_map_file="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr$_b37_to_b38_map.txt"
info_file="/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/impute_$_interval.snpstats"

library(dplyr)
library(data.table)

cis=fread(pergene)
cis2=cis%>%
  filter(qval<=0.05)

for(i in 1:22) {
  cat(paste0("\nReading in Chr", i))
  cat(paste0("\n", b37_b38_map_file))
  cat(paste0("\n", gsub(b37_b38_map_file, pattern="\\$", replacement=i)))
  cis_nom <- fread(data.table = F, paste0(prefix_nominal, i, ".csv"))
  #cis <- fread(data.table = F, paste0(prefix_pergene, i, ".csv"))
  map <- fread(data.table = F, paste0(prefix_bim, i, ".bim"))
  names(map) <- c("chr","variant_id","morg","pos_b38", "effect_allele", "other_allele")
  map <- select(map, -morg)
  b37_b38_map=fread(gsub(b37_b38_map_file, pattern="\\$", replacement=i))
  info=fread(gsub(info_file, pattern="\\$", replacement=i))

  b37_b38_map=b37_b38_map%>%
    mutate(rsid=ifelse(rsid==".", paste(chromosome, position, alleleA, alleleB, sep="_"), rsid))

  cis_nom2=cis_nom%>%
    left_join( map, by = "variant_id")%>%
    left_join(b37_b38_map%>%
      select(rsid.b38,rsid,position)%>%
      rename(variant_id_b37=rsid, pos_b37=position), by=c("variant_id"="rsid.b38"))%>%
    left_join(info%>%select(SNPID, information), by=c("variant_id_b37"="SNPID"))

  eSNPs=cis_nom2%>%
    left_join(cis2%>%select(Phenotype, pval_nominal_threshold), by=c("phenotype_id"="Phenotype"))%>%
    filter(pval_nominal<=pval_nominal_threshold)%>%
    select(-pval_nominal_threshold)
  
  numsigneSNPs=eSNPs%>%
    select(phenotype_id, variant_id)%>%
    distinct()%>%
    count(phenotype_id)%>%
    rename(num_sig_eSNPs=n)

  eGenes=cis2%>%
    left_join(map)%>%
    filter(chr==i)%>%
    left_join(numsigneSNPs, by=c("Phenotype"="phenotype_id"))
  
  fwrite(eGenes, file = paste0(prefix_results, i, "_significant_eGenes.csv"))
  fwrite(eSNPs, file = paste0(prefix_results, i, "_significant_eSNPs.csv"))
  fwrite(cis_nom2, file = paste0(prefix_nominal, i, ".annotated.csv"))

  rm(cis_nom, cis_nom2, map, eGenes, eSNPs)
}
              
              
