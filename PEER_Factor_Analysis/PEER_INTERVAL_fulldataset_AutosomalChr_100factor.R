# ----------------------------------------
# Running peer analysis 
# ----------------------------------------
# Gene expression data has 4143samples x 17964 genes
# Gene data TMM nornalised and each gene is inverse rank normalised
#subset the data below to only include autosomal chromosomes

#gene expression   -  matrix of gene expression data - columns are genes, rows are samples
#covariates matrix - columns are covariates are, e.g. Age, Sex, BMI -- Must not have missing data!
# Note the sample order must be the same for 'gene expression data' and 'covariates data', PEER won't match up the samples for you

# Covariates data includes 
# we including age, sex, BMI and blood cell traits: 
#BASO_PCT___RNA, EO_PCT___RNA, LYMPH_PCT___RNA, MONO_PCT___RNA,NEUT_PCT___RNA
#WBC_10_9_L___RNA, IRF_PCT___RNA, RET_PCT___RNA, HCT_PCT___RNA, HGB_g_dL___RNA, MCH_pg___RNA
#MCHC_g_dL___RNA, MCV_fL___RNA, RBC_10_12_L___RNA, RDW_SD_fL___RNA, MPV_fL___RNA, PCT_PCT___RNA, PDW_fL___RNA, PLT_F_10_9_L___RNA

#swap the swapped samples and remove samples without matching genotypes

library(peer)
setwd("~/06_AllSamples_Final_Analysis")
gene_expr <- read.csv("Input_Files/GeneExpr_PEER_TmmInvRankNormalised_swapsSwapped_mismatchRemoved.csv", row.names = 1, stringsAsFactors = FALSE)
covariates <-read.csv("Input_Files/Covariates_PEER_swapsSwapped_mismatchRemoved.csv", row.names = 1, stringsAsFactors = FALSE)
anno <- read.csv("Input_Files/GeneAnnotations.csv", stringsAsFactors = FALSE)

table(rownames(gene_expr) == rownames(covariates))

table(anno$gene_id == colnames(gene_expr))
anno1 <- anno[which(as.numeric(anno$chr) %in% c(1:22)),]

gene_expr <- gene_expr[,which(colnames(gene_expr) %in% anno1$gene_id)]  #19173 genes


gene_expr <- as.matrix(gene_expr)
covariates <- as.matrix(covariates)

test <- covariates[1:10, 1:10]
write.csv(test, "PEER/out_100Factor_AutosomalChr/testfile.csv")

# Run PEER

print("running PEER with genes from autosomal chromosomes - estimating 100 factors, run date 6th May, 2021")

# Set up model
model <- PEER()

# adds intercept term to PEER linear model fits

# adds intercept term to PEER linear model fits

PEER_setPhenoMean(model, gene_expr)
dim(PEER_getPhenoMean(model))

# adds additional factor to account for the mean expression
PEER_setAdd_mean(model, TRUE)

# PEER will find factors independent of these covariates
PEER_setCovariates(model, covariates)
dim(PEER_getCovariates(model))

# For automatic factor selection the PEER Nature Protocols paper recommends
# setting the number of factors to 25% of the sample size to a maximum of 100.
# We can then select the appropriate number of factors by looking at their variance.
nK <- min(75, round(nrow(gene_expr)/4))
nK
PEER_setNk(model, nK)
PEER_getNk(model)

# Increase to 10,000 if the next command gives a message about failing to converge
PEER_setNmax_iterations(model, 1000)

# Run PEER factor analysis
PEER_update(model)

# Get the hidden factors that PEER identified
factors <- PEER_getX(model)
dim(factors)
rownames(factors) <- rownames(gene_expr)
colnames(factors) <- c(colnames(covariates), "intercept", paste0("PEER", 1:nK))
saveRDS(factors, "PEER/out_100Factor_AutosomalChr/PEER_factor.rds")
write.table(factors,"PEER/out_100Factor_AutosomalChr/PEER_factors.txt",quote=F,sep="\t")

# Get the variance explained by each covariate and factor - we will use this to choose
# the number of factors to use in downstream analyses
factor_variance <- PEER_getAlpha(model)
factor_variance <- 1/factor_variance # Convert from inverse variance to variance
rownames(factor_variance) <- c(colnames(covariates), "intercept", paste0("PEER", 1:nK))
factor_variance <- factor_variance[-(1:(ncol(covariates)+1)),] # ignore covariates and intercept
saveRDS(factor_variance, "PEER/out_100Factor_AutosomalChr/PEER_factors_varaince_explained.rds")
write.table(factor_variance,"PEER/out_100Factor_AutosomalChr/PEER_factorvariance.txt",quote=F,sep="\t")

# Look at this plot: for downstream analysis you want to choose the first N factors
# before the variance drops to 0
pdf("PEER/out_100Factor_AutosomalChr/PEER_factor_kN_selection.pdf")
plot(factor_variance, type="l", ylab="Variance of PEER factor weights", xlab="", xaxt="n")
points(factor_variance, pch=19)
dev.off()

# Other information from the PEER analysis we don't need, but you might
# want to save anyway:
# The weights correspond to how much each factor/covariate affects each gene
weights <- PEER_getW(model)
colnames(weights) <- c(colnames(covariates), "intercept", paste0("PEER", 1:nK))
rownames(weights) <- colnames(gene_expr)
write.table(weights,"PEER/out_100Factor_AutosomalChr/PEER_weights.txt",quote=F,sep="\t")
