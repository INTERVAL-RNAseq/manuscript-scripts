#!/bin/bash

chr=$1

# module load ceuadmin/gcta/1.26.0
module load plink-1.9-gcc-5.4.0-sm3ojoi

# get summary statistics file for one gene
nominal_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cisNominal_chr${chr}.annotated.csv"
egene_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/tensorqtl_cis_MAF0.005_cis_chr${chr}_significant_eGenes.csv"
plink_prefix_path="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/genotypes/sorted/INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_chr${chr}_sorted"
list_output="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/cojooutput_tensorqtl_cis_MAF0.005_cisNominal_chr${chr}.outputs"
gcta="/home/ep620/rds/hpc-work/programs/gcta/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static"

genes=` awk 'BEGIN { FS = "," } ; NR>1{print $1}' ${egene_file} `

#gene="ENSG00000275329"

for gene in $genes
do
    cojoinput_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/cojoinput_tensorqtl_cis_MAF0.005_cisNominal_chr${chr}.annotated.${gene}.txt"
    output_file="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/cojooutput_tensorqtl_cis_MAF0.005_cisNominal_chr${chr}.annotated.${gene}.txt"
    list_snps="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/cojoinput_listsnps.${gene}.txt"
    cojoplink_prefix_path="/home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/INTERVAL_RNAseq_Phase1-3_imputed_b38_biallelic_MAF0.005_chr${chr}_sorted.${gene}"

    # columns should be SNP A1 A2 freq b se p N 
    variant=`awk -v gene=${gene} 'BEGIN { FS = "," } ;NR>1{if($1==gene) print $7}' $egene_file`
    awk -v gene=${gene} -v top=${variant} 'BEGIN { FS = "," } ;NR==1{print "SNP","A1","A2","freq","b","se","p","N"};NR>1{if($1==gene){if($2==top){print $2,$12,$13,$4,$8,$9,$7,"4732"}else{print $2,$12,$13,$4,$8,$9,"1","4732"}}}' $nominal_file > $cojoinput_file
    #awk -v gene=${gene} 'BEGIN { FS = "," } ;NR==1{print "SNP","A1","A2","freq","b","se","p","N"};NR>1{if($1==gene)print $2,$12,$13,$4,$8,$9,$7,"4732"}' $nominal_file > $cojoinput_file
    awk -v gene=${gene} 'BEGIN { FS = "," } ;NR>1{if($1==gene) print $2}' $nominal_file > $list_snps
    pval_threshold=`awk -v gene=${gene} 'BEGIN { FS = "," } ;NR>1{if($1==gene) print $18}' $egene_file`

    # select snps for the gene 
    plink \
        --bfile $plink_prefix_path \
        --extract ${list_snps} \
        --make-bed \
        --out ${cojoplink_prefix_path}

    # run cojo on this file
    $gcta  --bfile ${cojoplink_prefix_path} --chr ${chr} --maf 0.001 --cojo-file ${cojoinput_file} --cojo-slct --cojo-p ${pval_threshold}  --out ${output_file}
    #  --cojo-collinear 0.9   --cojo-actual-geno deprecated

    echo $gene

    rm ${cojoplink_prefix_path}.bed
    rm ${cojoplink_prefix_path}.fam
    rm ${cojoplink_prefix_path}.bim
    rm ${cojoplink_prefix_path}.nosex
    rm ${cojoplink_prefix_path}.log
    rm ${cojoinput_file}
    rm ${list_snps}
done


ls /home/ep620/rds/rds-jmmh2-projects/interval_rna_seq/analysis/06_AllSamples_Final_Analysis/results/cis_eQTLs/indep_112022/cojooutput_tensorqtl_cis_MAF0.005_cisNominal_chr${chr}.annotated.*.txt.jma.cojo > ${list_output}