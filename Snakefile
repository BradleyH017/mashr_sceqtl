########################################################################################################
################## Snakefile to run mashr on multiple conditions of eQTL / DE tests ####################
########################################################################################################
###################################### Author: Bradley Harris ##########################################
########################################## Date: 22/07/24 ##############################################
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

conditions = get_lines(config["conditions_file"])

# Load variables
outdir = config["outdir"] 
inputdir = config["input_dir"]
chrs = config["chrs"]
suffix = config["suffix"]
nchunks = config["nchunks"]

# define gather input function
def gather_nominal(wildcards):
    return expand(inputdir + "/{condition}/" + suffix + "{chr}.tsv",
                  condition=conditions,
                  chr=chrs)

# Define rule all
rule all: 
    output:
        outdir + "/mash_model/mash.rds"
        #outdir + "/input/sumstats_subset.rds"
        #expand(outdir + "/postmashr/{conditions}_chr{chr}_prediction.txt.gz", outdir = outdir, chr = chrs)

# Define execution rules:
rule aggregate_sumstats_filter:
    input:
        gather_nominal
    output:
        outdir + "/input/sumstats_subset.rds",
        expand(outdir + "/input/chunks/chunk{chunk}.rds", chunk = nchunks) # Not sure this will work 
    params:
        input_dir = config["input_dir"],
        conditions_file = config["conditions_file"],
        chrs = config["chrs"],
        suffix = config["suffix"],
        outdir = config["outdir"],
        nrand = config["nrand"], 
        fill_missing_beta = config["fill_missing_beta"],
        fill_missing_se = config["fill_missing_se"]
    singularity:
        "/software/hgi/softpack/installs/groups/otar2065//seurat5_v2/1-scripts/singularity.sif"
    resources:
    shell: 
        r"""
        mkdir -p {params.outdir}/input
        Rscript bin/001-gather_sumstats_all_nominal.r \
            --input_dir {params.input_dir} \
            --conditions_file {params.conditions_file} \ 
            --chrs {params.chrs} \
            --suffix {params.suffix} \
            --outdir {params.outdir} \
            --nrand {params.nrand} \
            --fill_missing_beta {params.fill_missing_se} \
            --fill_missing_se {params.fill_missing_se}
        """


rule gen_model:
    input:
        outdir + "/input/sumstats_subset.rds"
    output:
        outdir + "/mash_model/mash.rds"
    params:
        outdir = config["outdir"]
        lfsr = config["lfsr"]
    resources:
    singularity:
        "/software/hgi/softpack/installs/groups/otar2065//seurat5_v2/1-scripts/singularity.sif"
    shell:
        r"""
        mkdir -p {params.outdir}/mash_model
        Rscript 002-gen_model.r \
            --input
            --outdir {params.outdir}/mash_model \
            --lfsr {params.lfsr} \
            
        """

rule apply_model:
