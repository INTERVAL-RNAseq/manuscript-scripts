library(dplyr)
library(data.table)
library(biomaRt)
library(httr)
set_config(config(ssl_verifypeer = 0L))

listEnsemblArchives()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl", host="https://jan2020.archive.ensembl.org")
attributes = listAttributes(ensembl)
gene_ensembl=getBM(attributes=c('chromosome_name',
                    'ensembl_gene_id',
                    'external_gene_name',
                    'start_position', 
                    'end_position',
                    'strand',
                    'transcription_start_site',
                    'description',
                    'gene_biotype',
                    'hgnc_id',
                    'hgnc_symbol',
                    'entrezgene_id'
                    ),
mart = ensembl)%>%
    mutate(across(everything(), as.character))

gene_ensembl=gene_ensembl%>%
    filter(transcription_start_site==start_position | transcription_start_site==end_position)%>%
    mutate(hgnc_symbol=ifelse(hgnc_symbol=="", NA, hgnc_symbol))%>%
    mutate(hgnc_symbol=ifelse(hgnc_id=="", NA, hgnc_symbol))%>%
    rename(feature_id=ensembl_gene_id,
            chromosome=chromosome_name,
            start=start_position,
            end=end_position,
            feature_strand=strand,
            gene_name=external_gene_name, 
            TSS=transcription_start_site)%>%        
    dplyr::select(feature_id, chromosome, start, end, feature_strand, gene_name, TSS, gene_biotype, description, gene_biotype, hgnc_id, hgnc_symbol, entrezgene_id)

write.table(gene_ensembl, "~/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/annotation_file/Feature_Annotation_Ensembl_gene_ids_ensembl99.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# checking if all gene ids from raw rnaseq count matrix are annotated
initialdata=fread("~/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/Data/Final_QCd_Filtered/INTERVAL_FilteredSamplesGenes_swapsSwappedMismatchRemoved_TMMNorm_InverseRankedValues.csv")
# 19835 initialdata_genes
initialdata_genes=colnames(initialdata)[-1]
nomatch2=which(!initialdata_genes %in% gene_ensembl$ensembl_gene_id)
#0 which don't match
nomatch2=which(!initialdata_genes %in% gene_ensembl$feature_id)