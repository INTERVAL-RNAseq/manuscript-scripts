# Python script to run trans-eQTL mapping
import torch
import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

# Specify root directory for analysis
path = str(sys.argv[1]) 

# Specify output path (and prefix if desired). Can be set to any folder if you want it outside the analysis directory
outpath = path + "results/trans_eQTLs/"

# Get analysis chromosome from command line arguments, this in turn comes from SLURM array ID
# chr = str(sys.argv[2]) 

# Set up file paths
phenotype_bed_file=str(sys.argv[2])
covariates_file=str(sys.argv[3])
plink_prefix_path=str(sys.argv[4])
indepsignals_file=str(sys.argv[5])
trans_df_outfile=str(sys.argv[6])
print("output file is " + trans_df_outfile)

# phenotype_bed_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr17.bed.gz"
# covariates_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/covariates/INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis_affyID.txt"
# plink_prefix_path="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/genotypes/INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_chr17"
# indepsignals_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_summary/allresults_cojo.csv"

# Read in phenotypes
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)

# Read in covariates and make subset to only ids that are in the phenotype file
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
covariates_df = covariates_df[phenotype_df.columns].T

# Read in genotypes 
pr = genotypeio.PlinkReader(plink_prefix_path)

# load genotypes and variants into data frames
genotype_df = pd.DataFrame(pr.load_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# sorting variants with positions
variant_df2=variant_df.sort_values(["chrom", "pos"], ascending = (True, True))
variant_df=variant_df2
genotype_df2=genotype_df.reindex(variant_df.index)
genotype_df=genotype_df2

# Reading list of indep signals and selecting them
indepsignals=pd.read_csv(indepsignals_file)
#list_variants=list(set(indepsignals['variant_id']))
list_variants=list(set(indepsignals['SNP']))
genotype_df2=genotype_df.loc[genotype_df.index.isin(list_variants)]
variant_df2=variant_df.loc[variant_df.index.isin(list_variants)]

# Call trans-eQTLs
trans_df = trans.map_trans(genotype_df2, phenotype_df, covariates_df, return_sparse=True, pval_threshold=1, maf_threshold = 0.005, return_r2=True)
trans_df.to_csv(trans_df_outfile, index = False)

trans_df = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df2, window=5000000)
trans_df.to_csv(trans_df_outfile + ".filter_cis", index = False)

print("all done")
