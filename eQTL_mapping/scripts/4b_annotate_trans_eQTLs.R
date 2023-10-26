# Annotating trans-eQTL results

args=commandArgs(trailingOnly=TRUE)
trans_file=args[1]
chr=args[2]
prefix_bim=args[3]

library(data.table)
library(dplyr)

b37_b38_map_file="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr$_b37_to_b38_map.txt"
info_file="/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/impute_$_interval.snpstats"

b37_b38_map=fread(gsub(b37_b38_map_file, pattern="\\$", replacement=chr))
info=fread(gsub(info_file, pattern="\\$", replacement=chr))

b37_b38_map=b37_b38_map%>%
    mutate(rsid=ifelse(rsid==".", paste(chromosome, position, alleleA, alleleB, sep="_"), rsid))

map <- fread(data.table = F, paste0(prefix_bim, chr, ".bim"))
  names(map) <- c("chr","variant_id","morg","pos_b38", "effect_allele", "other_allele")
  map <- select(map, -morg)

trans=fread(trans_file)
trans2=trans%>%
    left_join( map, by = "variant_id")%>%
    left_join(b37_b38_map%>%
        select(rsid.b38,rsid,position)%>%
        rename(variant_id_b37=rsid, pos_b37=position), by=c("variant_id"="rsid.b38"))%>%
        left_join(info%>%select(position, information, A_allele, B_allele), by=c("pos_b37"="position", "effect_allele"="A_allele", "other_allele"="B_allele"))

trans2%>%filter(is.na(information))%>%head(1)

fwrite(trans2, paste(trans_file, ".annotated.csv", sep=""))