library("coloc")
library("stringr")
library("dplyr")
library("rlang")
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly = TRUE)
type = as.character(args[1]) # Phenotype code for COVID-HGI, i.e. "A2"
print(type)
prot = type

# Load map of genes to test (has eQTL within range of genome-wide significant SNP in COVID-19 HGI)
map = read.table(paste0("maps/map_",type,".txt"),
    header=FALSE,sep=" ")

# Case/control numbers for COVID-HGI
if (type=="A2") {
        cc=c(18152,1145546)
}
if (type=="B1") {
        cc=c(16512,71321)
}
if (type=="B2") {
        cc=c(44986,2356386)
}
if (type=="C2") {
        cc=c(159840,2782977)
}
cases_N=cc[1]
controls_N=cc[2]
total_N = cases_N+controls_N
case_PROP = cases_N/total_N

# Import COVID-HGI summary statistics
g_stats = read.table(paste0("gwas/COVID19_HGI_",type,"_ALL_leave_23andme_20220403.tsv.gz"),header=TRUE,sep="\t",comment.char = "")
colnames(g_stats)[1] = "CHR"

td = tempdir()
i = 1
l = list()
for (p in 1:nrow(map)) {
    chr = map$V2[p]
    splice = map$V1[p]
    print(paste0("Doing ",p,"/",nrow(map)))
    splice_p = map$V3[p]
    prot_p = 5e-8

    # Read in eQTL summary statistics for genes
    cmd = paste0("zcat /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq_n5188/eQTL_sumstats/tensorqtl_cis_MAF0.005_cisNominal_chr",
             chr,".annotated.csv.gz | ",'grep "^',splice,'"', " | awk -F',' ' {print $2,$12,$13,$4,$8,$9,$7,4732}' > ",td,"/sqtls")
    system(cmd)
    sqtls = read.table(paste0(td,"/sqtls"),header=FALSE,sep=" ")

    # Extract relevant SNPs
    write.table(unique(sqtls$V1),paste0(td,"/list_snps"),quote=FALSE,col.names=FALSE,row.names=FALSE)
    cmd = paste0("/software/team282/at24/conda/envs/plink/bin/plink ",
                     "--bfile /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq_n5188/alex/genotypes_sorted/out/INTERVAL_RNAseq_Phase1-3_imputed_filtered_b38_biallelic_MAF0.005_chr",chr," ",
                    "--extract ",td,"/list_snps ",
                    "--make-bed ",
                    "--out ",td,"/snps")
    system(cmd)

    # Calculate conditionally independent QTLs
    cmd = paste0("/software/team282/at24/gcta_v1.94.0Beta_linux_kernel_4_x86_64/gcta_v1.94.0Beta_linux_kernel_4_x86_64_static ",
                 "--bfile ",td,"/snps ",
                 "--cojo-file ",td,"/sqtls ",
                 "--cojo-slct --cojo-p ",splice_p," ",
                 "--out ",td,"/sqtls.cojo")
    system(cmd)
    pqtls = g_stats[which(g_stats$rsid %in% unique(sqtls$V1)),c("rsid","ALT","REF","all_meta_AF","all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p")]
    pqtls$N = total_N
    write.table(pqtls,paste0(td,"/pqtls"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep=" ")
    cmd = paste0("/software/team282/at24/gcta_v1.94.0Beta_linux_kernel_4_x86_64/gcta_v1.94.0Beta_linux_kernel_4_x86_64_static ",
                 "--bfile ",td,"/snps ",
                 "--cojo-file ",td,"/pqtls --diff-freq 0.2 ",
                 "--cojo-slct --cojo-p ",prot_p," ",
                 "--out ",td,"/pqtls.cojo")
    system(cmd,intern=TRUE)

    # Check these files exist, if not then no indep variants below significance threshold found
    if (file.exists(paste0(td,"/sqtls.cojo.jma.cojo"))) {
        cmd = paste0("tail -n+2 ",td,"/sqtls.cojo.jma.cojo | sort -gk 13 | awk '{print $2}' > ",td,"/sqtl_indep_snps.txt")
        system(cmd)
    } else {
        cat(paste0("EXITING: No indep sQTLs"))
        system(paste0("rm ",td,"/*"))
        next
    }
    if (file.exists(paste0(td,"/pqtls.cojo.jma.cojo"))) {
        cmd = paste0("tail -n+2 ",td,"/pqtls.cojo.jma.cojo | sort -gk 13 | awk '{print $2}' > ",td,"/pqtl_indep_snps.txt")
        system(cmd)
    } else {
        cat(paste0("EXITING: No indep pQTLs"))
        system(paste0("rm ",td,"/*"))
        next
    }

    # Make a matrix of eQTL and GWAS independent SNPs for pairwise coloc 
    splice_conds = cbind.data.frame(n=0,snp="unconditioned",file="sqtls")
    if (as.integer(unlist(str_split(system(paste0("wc -l ",td,"/sqtl_indep_snps.txt"),intern=TRUE)," "))[1])>1) {
        cmd = paste0('x=1; while read p; do cat ',td,'/sqtl_indep_snps.txt | grep -v "${p}" > ',td,'/snp_${x}_lead_${p}.txt; ',
                     '/software/team282/at24/gcta_v1.94.0Beta_linux_kernel_4_x86_64/gcta_v1.94.0Beta_linux_kernel_4_x86_64_static --bfile ',td,'/snps ',
                     '--cojo-file ',td,'/sqtls ',
                     '--cojo-cond ',td,'/snp_${x}_lead_${p}.txt ',
                     '--out ',td,'/INDEP_SPLICE_${x}_lead_${p}_out.txt; rm ',td,'/snp_${x}_lead_${p}.txt; x=$((x+1)); done < ',td,'/sqtl_indep_snps.txt')
        system(cmd)
        snames = dir(td,pattern = "INDEP_SPLICE_.*.cma.cojo")
        tmp = cbind.data.frame(as.data.frame(str_match(snames,"INDEP_SPLICE_([0-9]+)_lead_(.*)_out.txt.cma.cojo"))[,c(2,3)],file=dir(td,pattern = "INDEP_SPLICE_.*.cma.cojo"))
        colnames(tmp) = c("n","snp","file")
        rownames(tmp) = NULL
        splice_conds = rbind.data.frame(splice_conds,tmp)
    }
    splice_conds$file = paste0(td,"/",splice_conds$file)
    prot_conds = cbind.data.frame(n=0,snp="unconditioned",file="pqtls")
    if (as.integer(unlist(str_split(system(paste0("wc -l ",td,"/pqtl_indep_snps.txt"),intern=TRUE)," "))[1])>1) {
        cmd = paste0('x=1; while read p; do cat ',td,'/pqtl_indep_snps.txt | grep -v "${p}" > ',td,'/snp_${x}_lead_${p}.txt; ',
                     '/software/team282/at24/gcta_v1.94.0Beta_linux_kernel_4_x86_64/gcta_v1.94.0Beta_linux_kernel_4_x86_64_static --bfile ',td,'/snps ',
                     '--cojo-file ',td,'/pqtls ',
                     '--cojo-cond ',td,'/snp_${x}_lead_${p}.txt ',
                     '--out ',td,'/INDEP_PROT_${x}_lead_${p}_out.txt; rm ',td,'/snp_${x}_lead_${p}.txt; x=$((x+1)); done < ',td,'/pqtl_indep_snps.txt')
        system(cmd)
        pnames = dir(td,pattern = "INDEP_PROT_.*.cma.cojo")
        tmp = cbind.data.frame(as.data.frame(str_match(pnames,"INDEP_PROT_([0-9]+)_lead_(.*)_out.txt.cma.cojo"))[,c(2,3)],file=dir(td,pattern = "INDEP_PROT_.*.cma.cojo"))
        colnames(tmp) = c("n","snp","file")
        rownames(tmp) = NULL
        prot_conds = rbind.data.frame(prot_conds,tmp)
    }
    prot_conds$file = paste0(td,"/",prot_conds$file)

    # Perform pairwise coloc
    for (s in 1:nrow(splice_conds)) {
        if (s==1) {
            s_stats = read.table(splice_conds$file[s],header=FALSE,sep=" ")
            s_stats = cbind.data.frame(SNP=s_stats$V1,pos=NA,beta=s_stats$V5,vb=s_stats$V6)
            s_stats = s_stats[which(!duplicated(s_stats$SNP)),]
            rownames(s_stats) = s_stats$SNP
        } else {
            s_stats = read.table(splice_conds$file[s],header=TRUE,sep="\t")
            rownames(s_stats) = s_stats$SNP
            s_stats = s_stats[,c(2,3,11,12)]
            s_stats = s_stats[which(!is.na(s_stats$bC)),]
        }
        # Subset for if conditioned or indep// ie if s_i = 1 or >1
        for (p in 1:nrow(prot_conds)) {
            if (p==1) {
                p_stats = read.table(prot_conds$file[p],header=FALSE,sep=" ")
                p_stats = cbind.data.frame(SNP=p_stats$V1,pos=NA,beta=p_stats$V5,vb=p_stats$V6)
                p_stats = p_stats[which(!duplicated(p_stats$SNP)),]
                rownames(p_stats) = p_stats$SNP
            } else {
                p_stats = read.table(prot_conds$file[p],header=TRUE,sep="\t")
                rownames(p_stats) = p_stats$SNP
                p_stats = p_stats[,c(2,3,11,12)]
                p_stats = p_stats[which(!is.na(p_stats$bC)),]
            }
            in_both = intersect(s_stats$SNP,p_stats$SNP)
            p_stats = p_stats[in_both,]
            s_stats_tmp = s_stats[in_both,]
            # Subset for if conditioned or indep// ie if s_i = 1 or >1

            sqtls_coloc = list()
            sqtls_coloc$snp = as.character(s_stats_tmp[,1])
            if (s==1) {
                sqtls_coloc$position = NULL
            } else {
                sqtls_coloc$position = as.character(s_stats_tmp[,2])
            }
            sqtls_coloc$beta = as.double(s_stats_tmp[,3])
            sqtls_coloc$varbeta = as.double(s_stats_tmp[,4])^2
            sqtls_coloc$type = "quant"
            sqtls_coloc$sdY = 1
            sqtls_coloc$N = 4732
            check_dataset(sqtls_coloc)

            pqtls_coloc = list()
            pqtls_coloc$snp = as.character(p_stats[,1])
            if (p==1) {
                pqtls_coloc$position = NULL
            } else {
                pqtls_coloc$position = as.character(p_stats[,2])
            }
            pqtls_coloc$beta = as.double(p_stats[,3])
            pqtls_coloc$varbeta = as.double(p_stats[,4])^2
            pqtls_coloc$type = "cc"
            pqtls_coloc$s = case_PROP
            pqtls_coloc$N = total_N
            check_dataset(pqtls_coloc)

            res=coloc.abf(dataset1 = sqtls_coloc, dataset2 = pqtls_coloc)

            l[[i]] = cbind.data.frame(chrom=chr,prot=type,
                                      splice=splice,
                                      snp_1_n = splice_conds$n[s],
                                      snp_1_id = splice_conds$snp[s],
                                      snp_2_n = prot_conds$n[p],
                                      snp_2_id = prot_conds$snp[p],
                                      nsnps=as.character(res$summary['nsnps']),
                                      H0=as.double(res$summary['PP.H0.abf']),
                                      H1=as.double(res$summary['PP.H1.abf']),
                                      H2=as.double(res$summary['PP.H2.abf']),
                                      H3=as.double(res$summary['PP.H3.abf']),
                                      H4=as.double(res$summary['PP.H4.abf']))

            # If pass, move the sumstats for the splice_snp and prot_snp to /keep
            pass = ((l[[i]]$H4 + l[[i]]$H3) > 0.9) & ((l[[i]]$H4/l[[i]]$H3)>3)

            # Save independent summary statistics for portal if pass coloc threshold
            if (pass) {
                fn = paste0("keep/",splice,"__",splice_conds$n[1],"__",splice_conds$snp[1],".txt")
                if (!file.exists(fn)) {
                    to_save = read.table(splice_conds$file[1],header=FALSE,sep=" ")
                    to_save = subset(to_save,select= -c(2,3,4,8))
                    to_save$omic = "eqtl"
                    to_save$phen = splice
                    to_save$cond = splice_conds$snp[1]
                    #to_save$z_abs = abs(to_save$V5/to_save$V6)
                    write.table(to_save,fn,sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
                }

                fn = paste0("keep/",splice,"__",prot,"__",prot_conds$n[1],"__",prot_conds$snp[1],".txt")
                if (!file.exists(fn)) {
                    to_save = read.table(prot_conds$file[1],header=FALSE,sep=" ")
                    to_save = subset(to_save,select= -c(2,3,4,8))
                    to_save$omic = "covid"
                    to_save$phen = prot
                    to_save$cond = prot_conds$snp[1]
                    #to_save$z_abs = abs(to_save$V5/to_save$V6)
                    #to_save$V7 = 10^to_save$V7
                    write.table(to_save,fn,sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
                }

                if (s>1) {
                    fn = paste0("keep/",splice,"__",splice_conds$n[s],"__",splice_conds$snp[s],".txt")
                    if (!file.exists(fn)) {
                        to_save = read.table(splice_conds$file[s],header=TRUE,sep="\t")
                        to_save = to_save[,c("SNP","bC","bC_se","pC")]
                        to_save = to_save[which(!is.na(to_save$bC)),]
                        to_save$omic = "eqtl"
                        to_save$phen = splice
                        to_save$cond = splice_conds$snp[s]
                        #to_save$z_abs = abs(to_save$bC/to_save$bC_se)
                        write.table(to_save,fn,sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
                    }
                }

                if (p>1) {
                    fn = paste0("keep/",splice,"__",prot,"__",prot_conds$n[p],"__",prot_conds$snp[p],".txt")
                    if (!file.exists(fn)) {
                        to_save = read.table(prot_conds$file[p],header=TRUE,sep="\t")
                        to_save = to_save[,c("SNP","bC","bC_se","pC")]
                        to_save = to_save[which(!is.na(to_save$bC)),]
                        to_save$omic = "covid"
                        to_save$phen = prot
                        to_save$cond = prot_conds$snp[p]
                        #to_save$z_abs = abs(to_save$bC/to_save$bC_se)
                        write.table(to_save,fn,sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE)
                    }
                }
            }
                       
            i = i+1

        }
    }
    # Tidy up temp files
    system(paste0("rm ",td,"/*"))
}

# Collate results and save
ll = do.call(rbind.data.frame,l)
write.table(ll,paste0("out/",type,".tsv"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

