# Instructions for running mashr

## 1. Prep the input
Adjust and run bin/000-prep_input_from_TQTL_minimal.sh
This will save a list of file names for conditions to use in the analysis. NOTE: This is optimised for use with TensorQTL output

## 2. Adjust config
Adjust the config file (config.yaml) so that it points to the corret directory for dataset input. 

## 3. Run snakemake pipeline
Run Snakefile_001_prep_gene_model.smk. On LSF cluster this can be done by submitting submit_snakemake.sh
This will: 
    1. aggregate the results from different chromosomes (if optimised for using TensorQTL output)
    2. Convert the output to .h5 and .rds files for input to mashr. Also selecting the set of strong and random tests
    3. Generate the mashr model
    4. Apply this model to the test results of each condition and save

