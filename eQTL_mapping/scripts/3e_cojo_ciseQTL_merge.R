# Concatenating all conditional analysis results

library(data.table)
library(dplyr)

files=dir("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/", pattern=".jma.cojo")

all_results=data.frame()
count=0
for(file in files){
    data=fread(paste("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/", file, sep=""))
    data=data%>%
        mutate(filename=file)
    all_results=bind_rows(all_results, data)
    count=count+1
    print(count)
}

all_results=all_results%>%
    mutate(gene=gsub(filename, pattern="cojooutput_tensorqtl_cis_MAF0.005_cisNominal_chr..annotated.", replacement=""))%>%
    mutate(gene=gsub(gene, pattern="cojooutput_tensorqtl_cis_MAF0.005_cisNominal_chr...annotated.", replacement=""))%>%
    mutate(gene=gsub(gene, pattern=".txt.jma.cojo", replacement=""))

# add tss distance 
nominal_all=data.frame()
for(i in 1:22){
  nominal=fread(paste("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cisNominal_chr",i,".annotated.csv", sep=""))
  nominal_all=bind_rows(nominal_all, nominal)
}

all_results2=all_results%>%
  left_join(nominal_all%>%select(phenotype_id, variant_id, tss_distance), by=c("gene"="phenotype_id", "SNP"="variant_id"))
write.csv(all_results2,"/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_summary/allresults_cojo.csv", row.names=FALSE)