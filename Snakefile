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
input_dir = config_load["input_dir"]
print("input_dir:", input_dir)
suffix_q = config["suffix_q"]
print("suffix_q:", suffix_q)

# Define rule all
rule all: 
    input:
        #"results/input/unique_genes.txt"
        #"fastqtl_to_mash_output/merged_test_conditions.mash.rds"
        "results/output/model.rds"

# Define execution rules:
if config["qtl_method"] == "TensorQTL":
    rule aggregate_sumstats_per_condition_TQTL:
        input:
            input_dir + "/{condition}/" + suffix_q
        output:
            "results/input/merged_{condition}.tsv.gz",
            "results/input/tested_genes_{condition}.tsv"
        params:
            input_dir=config["input_dir"],
            suffix_nom=config["suffix_nom"],
            suffix_q=config["suffix_q"]
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
            echo {wildcards.condition}
            fnom1={params.input_dir}/{wildcards.condition}/{params.suffix_nom}1.tsv
            head $fnom1 -n 1 >> results/input/merged_{wildcards.condition}.tsv
            for c in {params.input_dir}/{wildcards.condition}/{params.suffix_nom}*.tsv; do
                tail -n +2 ${{c}} >> results/input/merged_{wildcards.condition}.tsv
            done

            # get the genes
            fq={params.input_dir}/{wildcards.condition}/{params.suffix_q}
            awk '{{print $1}}' $fq |  tail -n +2 >> results/input/tested_genes_{wildcards.condition}.tsv

            # compress the output
            gzip results/input/merged_{wildcards.condition}.tsv
            """

def gather_unique_genes_across_conditions(wildcards):
    return expand("results/input/tested_genes_{condition}.tsv",
                  condition=conditions)

rule get_unique_genes:
    input:
        gather_unique_genes_across_conditions
    output:
        "results/input/unique_tested_genes.tsv",
        "results/input/merged_test_conditions.list"
    params:
        conditions_file=config["conditions_file"]
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
        # Get unique genes
        temp_file=$(mktemp)
        find results/input -type f -name 'tested_genes*' -exec cat {{}} + > "$temp_file"
        sort "$temp_file" | uniq > results/input/unique_tested_genes.tsv
        rm "$temp_file"

        # Also make list file
        while read i; do
            echo merged_${{i}}.tsv.gz >> results/input/merged_test_conditions.list
        done <{params.conditions_file}
        """

def gather_merged_genes_across_conditions(wildcards):
    return expand("results/input/merged_{condition}.tsv.gz",
                  condition=conditions)

# TO DO: Add option to specify the beta, se and p-value columns from the input - in case it doesn't match
rule make_h5s:
    input:
        "results/input/unique_tested_genes.tsv",
        "results/input/merged_test_conditions.list",
        gather_merged_genes_across_conditions
    output:
        "fastqtl_to_mash_output/merged_test_conditions.mash.rds"
    params:
        function=config["function"],
        nrand=config["nrand"],
        strong_per_gene=config["strong_per_gene"]
    singularity:
        config["h5_singularity"]
    resources:
        mem=50000,
        queue='long',
        mem_mb=50000,
        mem_mib=50000,
        disk_mb=50000,
        tmpdir="tmp",
        threads=4
    shell:
        r"""
        sos run {params.function} \
            --data-list {input[1]} \
            --gene-list {input[0]} \
            -j 8 \
            --best-per-gene {params.strong_per_gene} \
            --random-snp-size {params.nrand}
        """

###### THEN run gen model??? May be able to have a list of options in config for the data-driven matrices types to use (inc. FLASH) ######

rule gen_model:
    input:
        "fastqtl_to_mash_output/merged_test_conditions.mash.rds"
    output:
        "results/output/model.rds"
    params:
        reference = config["reference"],
        include_matrices = config["include_matrices"],
    resources:
        mem=50000,
        queue='long',
        mem_mb=50000,
        mem_mib=50000,
        disk_mb=50000,
        tmpdir="tmp",
        threads=4
    singularity:
        "/software/hgi/softpack/installs/groups/otar2065//seurat5_v2/1-scripts/singularity.sif"
    shell:
        r"""
        # Make dir
        mkdir -p results/output/
        
        # Flatten the data driven matrices
        matrices_combined="$(echo "{params.include_matrices}" | sed 's/ /,/g')"
        echo $matrices_combined

        # Pass to script
        Rscript bin/001-gen_model.r -i {input} -o "results/output" -m $matrices_combined -r {params.reference}
        """

'''
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
'''
