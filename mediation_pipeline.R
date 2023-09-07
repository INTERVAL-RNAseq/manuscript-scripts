options(stringsAsFactors=FALSE)
library(mediation)
library(vcfR)
library(stringr)
library(diagram)
library(medflex)
p0=paste0

# Load sQTL summaries
sqd = read.table("../fdr_lt_5pc_allAnnots.tsv",sep="\t",header=TRUE,quote="")
rownames(sqd) = sqd$phenotype_id

# Load independent QTLs
all_results = read.table("../cojo_summary.tsv",sep="\t",header = TRUE)
all_results = all_results[which(all_results$gene %in% sqd$phenotype_id),]

# Load pairs to test from the coloc
load("../coloc/sqtl_out/coloc_pass_all.RData")
coloc_cs_pairs = coloc_cs[,c("splice","prot")]
coloc_cs_pairs = coloc_cs_pairs[which(!duplicated(coloc_cs_pairs)),]
coloc_cs_pairs$prot = str_replace_all(tolower(coloc_cs_pairs$prot),"\\.","")

# Create matrix of pairs to test
l = list()
i = 1
for (splice in unique(coloc_cs_pairs$splice)) {
    snps = unique(all_results[which(all_results$gene==splice),]$SNP)
    somas = unique(coloc_cs_pairs[which(coloc_cs_pairs$splice==splice),]$prot)
    chr = sqd[splice,]$splice_chr
    for (snp in snps) {
        for (soma in somas) {
            l[[i]] = data.frame(splice=splice,soma=soma,lead=snp,chr=chr)
            i = i+1
        }
    }
}
coloc_cs_pairs = do.call(rbind.data.frame,l)

# Load protein phenotypes
gwasqc = read.table("somalogic_qcgwas_4000.csv.gz",
                    sep=",",header=TRUE,row.names=1)


# Run per chromosome
args = commandArgs(trailingOnly=TRUE)
CHR = as.integer(args[1])
cat(paste0("Doing CHR: ",CHR,"\n"))
l = list()
l_N = 1

# Load sQTL residuals
s_expr_all=read.table(paste0("/residuals/",CHR,".tsv"),
              sep="\t",header=TRUE,quote="")

coloc_cs_pairs_chr = coloc_cs_pairs[which(coloc_cs_pairs$chr==CHR),]

# Load VCF for the independent lead SNPs
vcf <- read.vcfR(paste0("/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq_n5188/alex/genotypes_sorted/sqtl_lead_indep_snps/out/lead_sqtl_chr",CHR,".vcf"), verbose = FALSE )
vcf_trim = as.data.frame(vcf@gt[,c(2:ncol(vcf@gt))])
vcf_trim = vcf_trim[which(!duplicated(getFIX(vcf)[,"ID"])),]
rownames(vcf_trim) = getFIX(vcf)[,"ID"][which(!duplicated(getFIX(vcf)[,"ID"]))]

if (nrow(coloc_cs_pairs_chr)<1) {
    write.table("",paste0("medflex_indep_out/soma_sqtl_medflex_summary_chr",CHR,".tsv"),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
} else {
    for (N in 1:nrow(coloc_cs_pairs_chr)) {

            cat(paste0("\t- Doing ",N,"/",nrow(coloc_cs_pairs_chr),"\n"))

            phen = coloc_cs_pairs_chr$splice[N]
            variant = coloc_cs_pairs_chr$lead[N]
            som = coloc_cs_pairs_chr$soma[N]

            cat(paste0("\tPhen: ",phen,", Som: ",som,"\n"))

            s_expr = as.data.frame(t(s_expr_all[phen,vcf_include]))
            colnames(s_expr) = c("s_expr")

            all_expr = s_expr

            print(table(rownames(s_expr) %in% rownames(vcf_trim)))
            print(variant)
            all_expr$snp = vcf_trim[rownames(s_expr),c(variant),drop=FALSE]

            all_expr = all_expr[grep("/",all_expr$snp),]
            if (length(grep("0/0",all_expr$snp))>0) {
                all_expr[which(all_expr$snp=="0/0"),]$snp = 0
            }
            if (length(grep("0/1",all_expr$snp))>0) {
                all_expr[which(all_expr$snp=="0/1"),]$snp = 1
            }
            if (length(grep("1/1",all_expr$snp))>0) {
                all_expr[which(all_expr$snp=="1/1"),]$snp = 2
            }
            all_expr$snp = as.numeric(all_expr$snp)

            all_expr$p_expr = gwasqc[rownames(all_expr),som]
            all_expr = all_expr[which(!is.na(all_expr$p_expr)),]

            fit.totaleffect = lm(formula=p_expr~snp, data=all_expr)
            tot_b = fit.totaleffect$coefficients["snp"]
            tot_se = coef(summary(fit.totaleffect))[, "Std. Error"]["snp"]
            tot_p = coef(summary(fit.totaleffect))[,4]["snp"]

            fit.mediator = lm(formula=s_expr~snp, data=all_expr)
            medi_b = fit.mediator$coefficients["snp"]
            medi_se = coef(summary(fit.mediator))[, "Std. Error"]["snp"]
            medi_p = coef(summary(fit.mediator))[,4]["snp"]

            fit.dv=lm(formula=p_expr~snp+s_expr, data=all_expr)

            dv_snp_b = fit.dv$coefficients["snp"]
            dv_snp_se = coef(summary(fit.dv))[, "Std. Error"]["snp"]
            dv_snp_p = coef(summary(fit.dv))[,4]["snp"]

            dv_splice_b = fit.dv$coefficients["s_expr"]
            dv_splice_se = coef(summary(fit.dv))[, "Std. Error"]["s_expr"]
            dv_splice_p = coef(summary(fit.dv))[,4]["s_expr"]

            rownames(all_expr) = str_replace(rownames(all_expr),"INT_RNA","")

            expData=neImpute(p_expr~snp+s_expr,data=all_expr,family=gaussian)
            nMod=neModel(p_expr ~ snp0 + snp1,family=gaussian,expData=expData, se = "robust")

            cf=coef(summary(neEffdecomp(nMod)))
            ci=confint(neEffdecomp(nMod))

            object = neEffdecomp(nMod)
            est = c(nMod$neModelFit$coefficients[2],nMod$neModelFit$coefficients[3],nMod$neModelFit$coefficients[3]+nMod$neModelFit$coefficients[2])
            se <- sqrt(diag(vcov(object)))
            zvalue <- coef(object)/se
            pvalue <- 2 * pnorm(-abs(zvalue))

            ade_est = est[1]
            ade_se = se[1]
            ade_z = zvalue[1]#cf[1,3]
            ade_p = pvalue[1]#cf[1,4]
            ade_c1 = ci[1,1]
            ade_c2 = ci[1,2]

            acme_est = est[2]
            acme_se = se[2]
            acme_z = zvalue[2]#cf[2,3]
            acme_p = pvalue[2]#cf[2,4]
            acme_c1 = ci[2,1]
            acme_c2 = ci[2,2]

            te_est = est[3]
            te_se = se[3]
            te_z = zvalue[3]#cf[2,3]
            te_p = pvalue[3]#cf[2,4]
            te_c1 = ci[3,1]
            te_c2 = ci[3,2]

            pm_est = acme_est/te_est

            l[[l_N]] = cbind.data.frame(N=nrow(all_expr),chr=CHR,splice=phen,variant=variant,phen=som,
                             tot_b=tot_b,tot_se=tot_se,tot_p=tot_p,
                             medi_b=medi_b,medi_se=medi_se,medi_p=medi_p,
                             dv_snp_b=dv_snp_b,dv_snp_se=dv_snp_se,dv_snp_p=dv_snp_p,
                             dv_splice_b=dv_splice_b,dv_splice_se=dv_splice_se,dv_splice_p=dv_splice_p,
                             ade_c1,ade_c2,ade_est,ade_p,ade_se,ade_z,acme_c1,acme_c2,acme_est,acme_p,acme_se,acme_z,te_c1,te_c2,te_est,te_p,te_se,te_z,pm_est)


            l_N = l_N+1
    }

    # Collate output and save
    res = do.call(rbind.data.frame,l)
    write.table(res,paste0("medflex_indep_out/soma_sqtl_medflex_summary_chr",CHR,".tsv"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}
