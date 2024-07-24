#!/usr/bin/env bash
# Merging results across chromosomes from TensorQTL output
# bsub -o logs/gather_across_chroms-%J-%I-output.log -e logs/gather_across_chroms-%J-%I-error.log -q normal -G team152 -n 1 -M 9000 -a "memlimit=True" -R "select[mem>9000] rusage[mem=9000] span[hosts=1]" -J "gather_across_chroms[1-5]" < bin/005-run_SAIGE_1_2_3_chrom1.sh 

# Load condition from the jobID
repo=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/mashr_sceqtl # Only bit that should change
cd $repo
condition_list=${repo}/results/input/test_conditions.txt
params_file=bin/prep_scripts/run_params.txt
condition=$(head $condition_list -n ${LSB_JOBINDEX} | tail -n 1)
inputdir=$(head bin/prep_scripts/run_params.txt -n 1 | tail -n 1)
nchrs=$(head bin/prep_scripts/run_params.txt -n 2 | tail -n 1)
suffix_nom=$(head bin/prep_scripts/run_params.txt -n 3 | tail -n 1)
suffix_qval=$(head bin/prep_scripts/run_params.txt -n 4 | tail -n 1)
outdir=$(head bin/prep_scripts/run_params.txt -n 5 | tail -n 1)

# Aggregate
echo $condition
fnom1=${inputdir}/${condition}/${suffix_nom}1.tsv
head $fnom1 -n 1 >> ${outdir}/merged_${condition}.tsv
for c in $(seq 1 $nchrs); do
    echo $c
    fnom=${inputdir}/${condition}/${suffix_nom}${c}.tsv
    tail -n +2 "$fnom" >> ${outdir}/merged_${condition}.tsv
done

# get the genes
fq=${inputdir}/${condition}/${suffix_qval}
awk '{print $1}' $fq |  tail -n +2 >> ${outdir}/tested_genes_${condition}.tsv

# compress the output
gzip ${outdir}/merged_${condition}.tsv
