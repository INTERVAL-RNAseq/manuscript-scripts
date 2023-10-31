# Code repository for the INTERVAL RNA-seq manuscript

The RNA-sequencing data pre-processing was done with a Nextflow pipeline available [here](https://github.com/wtsi-hgi/nextflow-pipelines/blob/rna_seq_interval_5591/pipelines/rna_seq.nf).

| File/Folder   | Description |
| -------- | ------- |
| `QC_Processing_RNAseq_ReadCounts` | Code for the processing of RNA-seq read count data |
| `PEER_Factor_Analysis` | Code for PEER factor analysis |
| `eQTL_mapping` | Pipeline for the eQTL mapping |
| `splice_event_annotation.R` | Pipeline for the annotation of sQTLs |
| `independent_coloc_pipeline.R` | Example pipeline of pairwise independent colocaliztion of eQTLs with COVID-19 HGI summary statistics |
| `mediation_pipeline.R` | Example pipeline of mediation analysis of colocalized sQTL-pQTL (Somalogic) pairs |

