# eQTL mapping pipeline

The eQTL mapping pipeline was run with Snakemake ([here](https://snakemake.readthedocs.io/en/stable/) for more information).
R and R packages were installed locally. R packages needed for the pipeline are :
* dplyr
* data.table
* biomaRt
* httr

The pipeline contains :
* `eQTLmapping_snakefile` : Snakefile 
* `envs/`: folder that contains the environment file `envs/tensorQTL_env.yml` to run the eQTL mapping.
* `scripts/`: folder that contains scripts being used in the pipeline:
    * `3_0_make_annotation_file_updated_ensembl99.R`
    * `3_3b_pergene_pvalthresholds.py` 
    * `3_4a_map_trans_eQTLs_cisindep.py`
    * `3_1_make_tensorQTL_input_phase2.R`
    * `3_3c_call_cis_eSNPs.R`
    * `3_4c_annotate_trans_eQTLs.R`
    * `3_2_index_bed.sh`
    * `3_3d_cojo_ciseQTL.sh`
    * `3_3a_map_cis_eQTLs.py`
    * `3_3e_cojo_ciseQTL_merge.R`

Note: Some of the code was copied and adapted from Jonathan Marten [github page](https://github.com/JonMarten/RNAseq) from first stage analysis.
