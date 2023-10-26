# eQTL mapping pipeline

The eQTL mapping pipeline was run with Snakemake ([here](https://snakemake.readthedocs.io/en/stable/) for more information).

Modules for R, plink and tabix from the CSD3 cluster were loaded. 

R packages were installed locally. R packages needed for the pipeline are :
* dplyr
* data.table
* biomaRt
* httr

The pipeline contains :
* `eQTLmapping_snakefile` : Snakefile 
* `envs/`: folder that contains the environment file `envs/tensorQTL_env.yml` to run the eQTL mapping.
* `scripts/`: folder that contains scripts being used in the pipeline:
    * `0_make_annotation_file_updated_ensembl99.R`
    * `1_make_tensorQTL_input_phase2.R`
    * `2_index_bed.sh`
    * `3a_map_cis_eQTLs.py`
    * `3b_pergene_pvalthresholds.py` 
    * `3c_call_cis_eSNPs.R`
    * `3d_cojo_ciseQTL.sh`
    * `3e_cojo_ciseQTL_merge.R`
    * `4a_map_trans_eQTLs_cisindep.py`
    * `4c_annotate_trans_eQTLs.R`

    
Note: Some of the code was copied and adapted from Jonathan Marten [github page](https://github.com/JonMarten/RNAseq) from first stage analysis.
