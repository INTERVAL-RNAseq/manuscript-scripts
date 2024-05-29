options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
library("stringr")
library("dplyr")
library("rlang")
library("GenomicRanges")
library("GenomicFeatures")
library("rtracklayer")
})

# All we need to start is a table with one column (here $phenotype_id) of the sQTL phenotype IDs
sqd = read.table("INTERVAL_sQTL_ids.tsv", header=TRUE)
rownames(sqd) = sqd$phenotype_id

# Extract intron chromosome, strand, start, end, cluster ID
splice_info = as.data.frame(t(as.data.frame(strsplit(sqd$phenotype_id,":")))[,c(1,2,3,4)])
colnames(splice_info) = c("splice_chr","splice_start","splice_end","splice_clu")
splice_info$splice_strand = str_sub(sqd$phenotype_id, -1, -1)
sqd = cbind.data.frame(sqd,splice_info)

# Read in leafcutter exon input (i.e. Gencode 33) and convert to ranges
exons = read.table("out_id_v33_gtype.txt.gz",
                  sep="\t", header=TRUE)
exons$chr = substr(exons$chr,4,30)
exons$gene_id = substr(exons$gene_id,1,15)
exons_granges = as(paste0(exons$chr,":",exons$start,"-",exons$end,":",exons$strand),"GRanges")
genes = exons %>% group_by(gene_id) %>% summarise(start = min(start),end=max(end),chr=max(chr),strand=max(strand),gene_name=max(gene_name),gene_type=max(gene_type))
genes_granges = as(paste0(genes$chr,":",genes$start,"-",genes$end,":",genes$strand),"GRanges")
values(genes_granges) = DataFrame(gene_id=genes$gene_id,gene_name=genes$gene_name,gene_type=genes$gene_type)
genes = as.data.frame(genes)
rownames(genes) = genes$gene_id
genes_granges_us = as(paste0(genes$chr,":",genes$start,"-",genes$end,":","*"),"GRanges")
values(genes_granges_us) = DataFrame(gene_id=genes$gene_id,gene_name=genes$gene_name)

# Match splice to genes fully containing the excised region
start_ranges = as(paste0(sqd$splice_chr,":",sqd$splice_start,":*"),"GRanges")
ov_start=as.data.frame(findOverlaps(start_ranges,genes_granges_us))
sqd$start_match = FALSE
for (i in 1:nrow(sqd)) {
    ge = genes[ov_start[which(ov_start$queryHits==i),]$subjectHits,]$gene_id
    start_matches = c()
    for (g in ge) {
        start_matches=c(start_matches,g)
    }
    sqd$start_matches[i] = paste0(start_matches,collapse=",")
}
end_ranges = as(paste0(sqd$splice_chr,":",sqd$splice_end,":*"),"GRanges")
ov_end=as.data.frame(findOverlaps(end_ranges,genes_granges_us))
sqd$end_match = FALSE
for (i in 1:nrow(sqd)) {
    ge = genes[ov_end[which(ov_end$queryHits==i),]$subjectHits,]$gene_id
    end_matches=c()
    for (g in ge) {
        end_matches=c(end_matches,g)
    }
    sqd$end_matches[i] = paste0(end_matches,collapse=",")
}

sqd$keep = FALSE
sqd$keep_genes = ""
for (i in 1:nrow(sqd)) {
    # If both have a value
    start_matches = unlist(str_split(sqd$start_matches[i],","))
    end_matches = unlist(str_split(sqd$end_matches[i],","))
    sqd$keep[i] = any(start_matches %in% end_matches)
    sqd$keep_genes[i] = paste0(start_matches[which(start_matches %in% end_matches)],collapse=",")
}

# Merge above to create final gene annotation for splice events, defaulting to sense-strand gene if there's a hit
sqd$merge_id = ""
sqd$merge_name = ""
sqd$merge_type = ""
sqd$merge_OS = FALSE
# Only iterate over rows in which keep_genes is not blank, therefore informative.
non_blank <- which(sqd$keep_genes != "")
for (i in non_blank) {
    mergeids = unlist(str_split(sqd[i,]$keep_genes,","))
    ## Order mergeids and supp by protein coding first, then others? (as other usually lncRNA)
    types = genes[mergeids,]$gene_type
    mergeids = mergeids[order(types=="protein_coding",decreasing = TRUE)]
    types = types[order(types=="protein_coding",decreasing = TRUE)]
    samestrand = c()
    for (m in mergeids) {
        samestrand = c(samestrand,genes[m,]$strand==sqd[i,]$splice_strand)
    }
    if (sum(samestrand)>0) {
        sqd$merge_id[i] = paste0(mergeids[samestrand],collapse=",")
        sqd$merge_name[i] = paste0(genes[mergeids[samestrand],]$gene_name,collapse=",")
        sqd$merge_type[i] = paste0(genes[mergeids[samestrand],]$gene_type,collapse=",")
    } else {
        sqd$merge_id[i] = paste0(mergeids," (OS)",collapse=",")
        sqd$merge_name[i] = paste0(genes[mergeids,]$gene_name, " (OS)",collapse=",")
        sqd$merge_type[i] = paste0(genes[mergeids,]$gene_type,collapse=",")
        sqd$merge_OS[i] = TRUE
    }
}

# Load full Gencode V33 GTF
gen33 <- rtracklayer::import('gencode.v33.annotation.gtf.gz')
gen33=as.data.frame(gen33)

sqd$range = paste0(sqd$splice_chr,":",sqd$splice_start,"-",sqd$splice_end,":",sqd$splice_strand)
# Whether spliced region/intron overlaps an exon
sqd$ov_exon = countOverlaps(as(sqd$range,"GRanges"),exons_granges,minoverlap=2,maxgap=-1)>0

# Whether the intron boundaries match known exons
exons$chrstart=""
exons$chrend=""
exons[which(exons$strand=="+"),]$chrstart = paste0(exons[which(exons$strand=="+"),]$chr,":",exons[which(exons$strand=="+"),]$start) # Different for neg strand?
exons[which(exons$strand=="+"),]$chrend = paste0(exons[which(exons$strand=="+"),]$chr,":",exons[which(exons$strand=="+"),]$end)
exons[which(exons$strand=="-"),]$chrstart = paste0(exons[which(exons$strand=="-"),]$chr,":",exons[which(exons$strand=="-"),]$end) # Different for neg strand?
exons[which(exons$strand=="-"),]$chrend = paste0(exons[which(exons$strand=="-"),]$chr,":",exons[which(exons$strand=="-"),]$start)
sqd$p5_exon = paste0(sqd$splice_chr,":",sqd$splice_start) %in% exons$chrend
sqd$p3_exon = paste0(sqd$splice_chr,":",sqd$splice_end) %in% exons$chrstart
sqd$p5_p3 = sqd$p5_exon & sqd$p3_exon
sqd$p5_or_p3 = sqd$p5_exon | sqd$p3_exon

# Find overlap of spliced region/introns with coding exons 
cds = gen33[which(gen33$gene_type=="protein_coding" & gen33$type=="CDS" & gen33$transcript_type=="protein_coding"),
          c("seqnames","start","end","strand","gene_name","gene_id")]
colnames(cds)[1] = "chr"
cds = cds[which(cds$chr %in% paste0("chr",1:22)),]
cds$chr = str_replace(cds$chr,"chr","")
cds$gene_id = substr(cds$gene_id,1,15)
cds_granges = as(paste0(cds$chr,":",cds$start,"-",cds$end,":",cds$strand),"GRanges")
sqd$ov_cds = countOverlaps(as(sqd$range,"GRanges"),cds_granges,minoverlap=2,maxgap=-1)>0

# Import and overlap UniProt/Pfam domain annotations (downloaded as .bed from UCSC genome browser)
valid_chr = c("chrX","chrY",paste0("chr",1:22))
unipLocTransMemb = read.table('ucsc_10feb22/unipLocTransMemb.bed', sep="\t", quote="")
unipLocTransMemb = unipLocTransMemb[which(unipLocTransMemb$V1 %in% valid_chr),]
unipLocTransMemb = as(paste0(unipLocTransMemb$V1,":",unipLocTransMemb$V2,"-",unipLocTransMemb$V3,":",unipLocTransMemb$V6),"GRanges")
unipLocTransMemb = reduce(unipLocTransMemb)
values(unipLocTransMemb) = DataFrame(type="Transmembrane")
unipLocCytopl = read.table('ucsc_10feb22/unipLocCytopl.bed', sep="\t", quote="")
unipLocCytopl = unipLocCytopl[which(unipLocCytopl$V1 %in% valid_chr),]
unipLocCytopl = as(paste0(unipLocCytopl$V1,":",unipLocCytopl$V2,"-",unipLocCytopl$V3,":",unipLocCytopl$V6),"GRanges")
unipLocCytopl = reduce(unipLocCytopl)
values(unipLocCytopl) = DataFrame(type="Cytoplasm")
unipLocExtra = read.table('ucsc_10feb22/unipLocExtra.bed', sep="\t", quote="")
unipLocExtra = unipLocExtra[which(unipLocExtra$V1 %in% valid_chr),]
unipLocExtra = as(paste0(unipLocExtra$V1,":",unipLocExtra$V2,"-",unipLocExtra$V3,":",unipLocExtra$V6),"GRanges")
unipLocExtra = reduce(unipLocExtra)
values(unipLocExtra) = DataFrame(type="Extracellular")
unipLocSignal = read.table('ucsc_10feb22/unipLocSignal.bed', sep="\t", quote="")
unipLocSignal = unipLocSignal[which(unipLocSignal$V1 %in% valid_chr),]
unipLocSignal = as(paste0(unipLocSignal$V1,":",unipLocSignal$V2,"-",unipLocSignal$V3,":",unipLocSignal$V6),"GRanges")
unipLocSignal = reduce(unipLocSignal)
values(unipLocSignal) = DataFrame(type="SignalP")
sqd$unipLocTransMemb = countOverlaps(as(paste0("chr",sqd$range),"GRanges"),unipLocTransMemb,minoverlap=2,maxgap=-1)>0
sqd[which(!sqd$ov_cds),]$unipLocTransMemb=FALSE
sqd$unipLocCytopl = countOverlaps(as(paste0("chr",sqd$range),"GRanges"),unipLocCytopl,minoverlap=2,maxgap=-1)>0
sqd[which(!sqd$ov_cds),]$unipLocCytopl=FALSE
sqd$unipLocExtra = countOverlaps(as(paste0("chr",sqd$range),"GRanges"),unipLocExtra,minoverlap=2,maxgap=-1)>0
sqd[which(!sqd$ov_cds),]$unipLocExtra=FALSE
sqd$unipLocSignal = countOverlaps(as(paste0("chr",sqd$range),"GRanges"),unipLocSignal,minoverlap=2,maxgap=-1)>0
sqd[which(!sqd$ov_cds),]$unipLocSignal=FALSE

unipDomainT = read.table('ucsc_10feb22/unipDomain.bed', sep="\t", quote="")
unipDomainT = unipDomainT[which(unipDomainT$V1 %in% valid_chr),]
unipDomain = as(paste0(unipDomainT$V1,":",unipDomainT$V2,"-",unipDomainT$V3,":",unipDomainT$V6),"GRanges")
values(unipDomain) = DataFrame(type=unipDomainT$V27)
ucscGenePfamT = read.table('ucsc_10feb22/ucscGenePfam.bed', sep="\t", quote="")
ucscGenePfamT = ucscGenePfamT[which(ucscGenePfamT$V1 %in% valid_chr),]
ucscGenePfam = as(paste0(ucscGenePfamT$V1,":",ucscGenePfamT$V2,"-",ucscGenePfamT$V3,":",ucscGenePfamT$V6),"GRanges")
values(ucscGenePfam) = DataFrame(type=ucscGenePfamT$V4)

x = as(paste0("chr",sqd$range),"GRanges")
hits=as.data.frame(findOverlaps(x,unipDomain,minoverlap=2,maxgap=-1))
doms = c()
for (i in 1:length(sqd$range)) {
    sub = hits[which(hits$queryHits==i),]
    if (dim(sub)[1]>0) {
        doms = c(doms,paste0(sort(unique(unipDomain$type[sub$subjectHits])),collapse = ","))
    } else {
        doms = c(doms,"-")
    }
}
sqd$unipDomain = doms
sqd[which(!sqd$ov_cds),]$unipDomain="-"
                                  
x = as(paste0("chr",sqd$range),"GRanges")
hits=as.data.frame(findOverlaps(x,ucscGenePfam,minoverlap=2,maxgap=-1))
doms = c()
for (i in 1:length(sqd$range)) {
    sub = hits[which(hits$queryHits==i),]
    if (dim(sub)[1]>0) {
        doms = c(doms,paste0(sort(unique(ucscGenePfam$type[sub$subjectHits])),collapse = ","))
    } else {
        doms = c(doms,"-")
    }
}
sqd$pfamDomain = doms
sqd[which(!sqd$ov_cds),]$pfamDomain="-"

sqd = sqd[,c('phenotype_id','splice_chr','splice_start','splice_end','splice_clu','splice_strand',
             'keep','merge_id','merge_name','merge_type','merge_OS','range','ov_exon',
             'p5_exon','p3_exon','p5_p3','p5_or_p3','ov_cds',
             'unipLocTransMemb','unipLocCytopl','unipLocExtra','unipLocSignal','unipDomain','pfamDomain')]


write.table(sqd,"sqtl_phenotypes_annotated.tsv",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

