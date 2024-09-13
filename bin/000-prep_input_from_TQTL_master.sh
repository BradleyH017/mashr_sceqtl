# Generating input conditions from TensorQTL output
# Following: https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb (need to also grab this file)
# https://stephenslab.github.io/gtexresults/fastqtl2mash.html
repo=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/mashr_sceqtl
inputdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/ashg_eqtl/rectum_label/results/TensorQTL_eQTLS
outdir=${repo}/results/input
nchrs=4
suffix_qval=OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv
suffix_nom=OPTIM_pcs/base_output__base/cis_nominal1.cis_qtl_pairs.chr
mkdir -p $outdir

# Define conditions to include in analysis
rm ${outdir}/test_conditions.txt
for dir in ${inputdir}/*; do
    f=${dir}/OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv
    if [ -f $f ]; then
        basedir=$(basename "$dir")
        echo "Adding: " ${basedir}
        echo ${basedir} >> ${outdir}/test_conditions.txt
    fi
done

# group results across chromsomes into a single file, gzip this and extract the list of tested genes
cd $repo
nconditions=$(wc -l < ${outdir}/test_conditions.txt)
# Echo params into a file that can be read by the script
rm bin/prep_scripts/run_params.txt # Remove if already present
echo $inputdir >> bin/prep_scripts/run_params.txt
echo $nchrs >> bin/prep_scripts/run_params.txt
echo $suffix_nom >> bin/prep_scripts/run_params.txt
echo $suffix_qval >> bin/prep_scripts/run_params.txt
echo $outdir >> bin/prep_scripts/run_params.txt

# Run a script to aggregate these files in parallel
mkdir -p logs
bsub -o logs/gather_across_chroms-%J-%I-output.log -e logs/gather_across_chroms-%J-%I-error.log -q normal -G team152 -n 1 -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -J "gather_across_chroms[1-${nconditions}]" <  bin/prep_scripts/000-001_merge_across_chromosomes_TQTL.sh

# Get a list of uniq genes
cd $repo
temp_file=$(mktemp)
find "$outdir" -type f -name 'tested_genes*' -exec cat {} + > "$temp_file"
sort "$temp_file" | uniq > ${outdir}/unique_tested_genes.tsv
rm "$temp_file"

# Add path to tested conditions to that the merged files are findable
while read i; do
    echo merged_${i}.tsv.gz >> ${outdir}/merged_test_conditions.list
done <${outdir}/test_conditions.txt

# Now run the conversion and prep of mash input files
# beta, se, pval = 8,9,7 (same as default)

# Installing the docker
module load ISG/singularity/3.11.4 
#export SINGULARITY_CACHEDIR=${repo}/singularity_cache
#mkdir -p ${repo}/singularity_cache
#singularity pull docker://gaow/hdf5tools

# Running
singularity shell -B /lustre -B /software $repo/hdf5tools_latest.sif
# Make sure we are where we want the output files to be
cd /lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/mashr_sceqtl/results/input_4cond_chr1
func=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/mashr_sceqtl/bin/fastqtl_to_mash.ipynb
# --best-per-gene 1 = Extracts the top eQTL per gene for 'strong' set
# --random-snp-size = Selects the number of random SNPs for the analysis
sos run $func \
    --data-list merged_test_conditions.list \
    --gene-list unique_tested_genes.tsv \
    -j 8 \
    --best-per-gene 0 \
    --random-per-gene -1 \
    --random-snp-size 200000