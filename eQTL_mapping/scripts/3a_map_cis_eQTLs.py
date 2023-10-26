# Python script to run cis-eQTL mapping
import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import gc

# Specify root directory for analysis
path = str(sys.argv[1]) 

# Specify output path (and prefix if desired). Can be set to any folder if you want it outside the analysis directory
outpath = path + "results/cis_eQTLs/"

# Get analysis chromosome from command line arguments, this in turn comes from SLURM array ID
chr = str(sys.argv[2]) 

# Set up file paths
phenotype_bed_file=str(sys.argv[3])
covariates_file=str(sys.argv[4])
plink_prefix_path=str(sys.argv[5])   

#cis_df_outfile=str(sys.argv[6]) 
cisnom_df_outprefix=str(sys.argv[6]) 

# chr=1
# phenotype_bed_file="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr1.bed.gz"
# covariates_file="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/covariates/INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis_affyID.txt"
# plink_prefix_path="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/genotypes/INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_chr1"
# cisnom_df_outprefix="/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/temp_chr1"

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

del variant_df2
del genotype_df2
gc.collect()

# Cis nominal mapping
cisnom_df = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, cisnom_df_outprefix, covariates_df, window=1000000, maf_threshold=0.005)
print("cis nominal mapping done")
cisnom_df2 = pd.read_parquet(cisnom_df_outprefix + ".cis_qtl_pairs." + chr + ".parquet")
cisnom_df2.to_csv(cisnom_df_outprefix + ".csv", index = False)


