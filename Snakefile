########################################################################################################
# Description: Snakefile to run mashr on multiple conditions of eQTL tests
# Author: Bradley Harris 
# Date: 22/07/24 
########################################################################################################

# load config
configfile: "config.yaml"

# Load the input conditions. These are to be used as wildcards
import yaml

configfile_yaml = open("config.yaml").read()
config_load = yaml.load(configfile_yaml, Loader=yaml.FullLoader)

def get_lines(conditions_file):
    with open(conditions_file, 'r') as f:
        values = [line.strip() for line in f]
    return values

conditions = get_lines(config_load["conditions_file"])
print("Conditions:", conditions)
inputdir = config_load["input_dir"]
print("input_dir:", input_dir)
suffix_q = config["suffix_q"]
print("suffix_q:", suffix_q)



# Define rule all
rule all: 
    input:
        expand(outdir + "/prep/{condition}/merged_nominal.tsv", condition = config["condition"])
        #outdir + "/input/sumstats_subset.rds"
        #expand(outdir + "/postmashr/{conditions}_chr{chr}_prediction.txt.gz", outdir = outdir, chr = chrs)


# Define execution rules:
rule aggregate_sumstats_per_condition:
    input:
        input_dir + "/{condition}/" + suffix_q
    output:
        "results/input/merged_{condition}.tsv.gz",
        "results/input/tested_genes_{condition}.tsv"
    params:
        input_dir=config["input_dir"],
        suffix_nom=config["suffix_nom"]
    resources:
        mem=5000,
        queue='normal',
        mem_mb=5000,
        mem_mib=5000,
        disk_mb=5000,
        tmpdir="tmp",
        threads=4
    shell: 
        r"""
        # Get the nominal
        echo ${wildcards.condition}
        fnom1={params.input_dir}/${condition}/{params.suffix_nom}1.tsv
        head $fnom1 -n 1 >> results/input/merged_${wildcards.condition}.tsv
        for c in {params.input_dir}/${condition}/{params.suffix_nom}*.tsv; do
            tail -n +2 ${{c}} >> results/input/merged_${wildcards.condition}.tsv
        done

        # get the genes
        fq={input]}
        awk '{print $1}' $fq |  tail -n +2 >> results/input/tested_genes_${params.condition}.tsv

        # compress the output
        gzip results/input/merged_${wildcards.condition}.tsv
        """

###### WRITE RULE TO DO THE CALCULATION OF UNIQUE GENES, AND ADDITION OF THE .LIST file for the input of these conditions

###### WRITE SCRIPT TO RUN THE WEIRD AGGREGATION - AND USE THE CUSTOM ENVIRONMENT TO EXTRACT RANDOM AND STRONG ######
###### Need to add options for number of random genes and --best-per-gene (1 = best eqtl only and is most appropriate) ######
###### THEN run gen model??? May be able to have a list of options in config for the data-driven matrices types to use (inc. FLASH) ######

rule gen_model:
    input:
        outdir + "/input/sumstats_subset.rds"
    output:
        outdir + "/mash_model/mash.rds"
    params:
        outdir = config["outdir"],
        reference = config["reference"],
        lfsr = config["lfsr"]
    resources:
    singularity:
        "/software/hgi/softpack/installs/groups/otar2065//seurat5_v2/1-scripts/singularity.sif"
    shell:
        r"""
        mkdir -p {params.outdir}/mash_model
        Rscript 002-gen_model.r \
            --input {input[0]} \
            --outdir {params.outdir}/mash_model \
            --reference {params.reference} \
            --lfsr {params.lfsr} 
        """

rule apply_model:
    input:
        expand(outdir + "/input/chunks/chunk{chunk}.rds", chunk = nchunks),
        outdir + "/mash_model/mash.rds"
    output:
        expand(outdir + "/output/chunks/chunk{chunk}.rds", chunk = nchunks)
    params:
        outdir = config["outdir"]
    resources:

    singularity:
        "/software/hgi/softpack/installs/groups/otar2065//seurat5_v2/1-scripts/singularity.sif"
    shell:
        r"""
        mkdir -p {params.outdir}/output/chunks
        Rscript
        """

