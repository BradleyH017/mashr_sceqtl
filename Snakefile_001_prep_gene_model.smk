########################################################################################################
# Description: Snakefile to prep input for mashr, then generate the model
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
input_dir = config_load["input_dir"]
suffix_q = config["suffix_q"]

# Define rule all
rule all: 
    input:
        #"results/input/unique_genes.txt"
        #"fastqtl_to_mash_output/merged_test_conditions.mash.rds"
        #"results/input/fastqtl_to_mash_output/merged_test_conditions.mash.rds","results/input/all_qval.tsv"
        "results/output/model.rds"

# Define execution rules:
if config["qtl_method"] == "TensorQTL":
    rule aggregate_sumstats_per_condition_TQTL:
        input:
            input_dir + "/{condition}/" + suffix_q
        output:
            "results/input/merged_{condition}.tsv.gz",
            "results/input/tested_genes_{condition}.tsv",
            "results/input/merged_{condition}_qval.tsv"
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

            # Copy the qval files here too, with new column for 
            cp $fq results/input/merged_{wildcards.condition}_qval.tsv
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

def gather_qval_across_conditions(wildcards):
    return expand("results/input/merged_{condition}_qval.tsv",
                  condition=conditions)

rule make_single_qval_file:
    input:
        gather_qval_across_conditions
    output:
        "results/input/all_qval.tsv"
    params:
        input_dir=config["input_dir"],
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
        temp_file=$(mktemp)
        for f in results/input/merged_*_qval.tsv; do
            echo $f
            if [ ! -f {output} ]; then
                echo "Getting header"
                head $f -n 1 >> $temp_file
                colname=condition
                awk -v value="$colname" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$temp_file" >> {output}
            fi
            base=$(basename "$f")
            awk -v value="$base" 'BEGIN{{FS=OFS="\t"}} {{print $0, value}}' "$f" | tail -n +2 >> {output}
        done
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
        "results/input/fastqtl_to_mash_output/merged_test_conditions.mash.rds",
        "results/input/fastqtl_to_mash_output/merged_test_conditions.h5"
    params:
        function=config["function"],
        random_per_gene=config["random_per_gene"],
        nrand=config["nrand"],
        strong_per_gene=config["strong_per_gene"]
    singularity:
        config["h5_singularity"]
    resources:
        mem=5000,
        queue='long',
        mem_mb=5000,
        mem_mib=5000,
        disk_mb=5000,
        tmpdir="tmp",
        threads=4
    shell:
        r"""
        sos run {params.function} \
            --data-list {input[1]} \
            --gene-list {input[0]} \
            -j 8 \
            --best-per-gene {params.strong_per_gene} \
            --random-per-gene {params.random_per_gene} \
            --random-snp-size {params.nrand} \
            --cwd results/input/fastqtl_to_mash_output
        
        #mv fastqtl_to_mash_output results/input/
        """

rule gen_model:
    input:
        "results/input/fastqtl_to_mash_output/merged_test_conditions.h5",
        "results/input/fastqtl_to_mash_output/merged_test_conditions.mash.rds",
        "results/input/all_qval.tsv"
    output:
        "results/output/model.rds",
        "results/output/complete_flag.txt"
    params:
        max_missing = config["max_missing"],
        replace_beta = config["replace_beta"],
        replace_se = config["replace_se"],
        lfsr_strong = config["lfsr_strong"],
        use_complete_random = config["use_complete_random"],
        nrand = config["nrand"],
        reference = config["reference"],
        include_matrices = config["include_matrices"],
        dd_matrix_version = config["dd_matrix_version"]
    resources:
        mem=50000,
        queue='long',
        mem_mb=50000,
        mem_mib=50000,
        disk_mb=50000,
        tmpdir="tmp",
        threads=4
    conda:
        "sm_mashr"
    shell:
        r"""
        # Make dir
        mkdir -p results/output/
        mkdir -p results/output/model_summary/predicted_results

        
        # Flatten the data driven matrices
        matrices_combined="$(echo "{params.include_matrices}" | sed 's/ /,/g')"
        echo $matrices_combined

        # Pass to script
        Rscript bin/001-gen_model_merged_h5.r -ih {input[0]} -im {input[1]} -mm {params.max_missing} -rb {params.replace_beta} -rs {params.replace_se} -qv {input[2]} -ls {params.lfsr_strong} -ucr {params.use_complete_random} -nr {params.nrand} -o "results/output" -m $matrices_combined -r {params.reference} -d {params.dd_matrix_version}
        """
