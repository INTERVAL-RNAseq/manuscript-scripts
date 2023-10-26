# Python script to run cis-eQTL mapping
import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis

# Set up file paths
phenotype_bed_file=str(sys.argv[1])
covariates_file=str(sys.argv[2])
plink_prefix_path=str(sys.argv[3])   
cis_df_outfile=str(sys.argv[4]) 

# phenotype_bed_file="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr"
# covariates_file="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/covariates/INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis_affyID.txt"
# plink_prefix_path="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/genotypes/INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_chr"
# cis_df_outfile="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cisPerGene_ALLchr.csv"

autosomal_chromosomes=[str(i) for i in range(1,23)]
cis_total_df=[]
for chrom in autosomal_chromosomes:
# Read in phenotypes
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file + chrom + ".bed.gz")
    # Read in covariates and make subset to only ids that are in the phenotype file
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    covariates_df = covariates_df[phenotype_df.columns].T
    # Read in genotypes 
    pr = genotypeio.PlinkReader(plink_prefix_path + chrom)
    # load genotypes and variants into data frames
    genotype_df = pd.DataFrame(pr.load_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    # sorting variants with positions
    variant_df2=variant_df.sort_values(["chrom", "pos"], ascending = (True, True))
    variant_df=variant_df2
    genotype_df2=genotype_df.reindex(variant_df.index)
    genotype_df=genotype_df2
    # Cis gene-level mapping
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, window=1000000, maf_threshold=0.005)
    cis_total_df.append(cis_df)

cis_total_df = pd.concat(cis_total_df)

print("gene level mapping done")
tensorqtl.calculate_qvalues(cis_total_df, fdr=0.05, qvalue_lambda=0) # Lambda of 0 is equivalent to BH correction. Prevents crashes for chr 9, 18, 22. 
print("q values computation done")
cis_total_df.to_csv(cis_df_outfile, index=True, index_label = "Phenotype")
