# Generating input conditions from TensorQTL output
# Following: https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb (need to also grab this file)
# https://stephenslab.github.io/gtexresults/fastqtl2mash.html
repo=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/mashr_sceqtl
inputdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_eqtl/2024_07_07-freeze005_TI_base_results/TensorQTL_eQTLS
outdir=${repo}/results/input
suffix_qval=OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv
suffix_nom=OPTIM_pcs/base_output__base/cis_nominal1.cis_qtl_pairs.chr
mkdir -p $outdir

# Define conditions to include in analysis
mkdir -p ${repo}/results/input
rm ${outdir}/test_conditions.txt
for dir in ${inputdir}/*; do
    f=${dir}/OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv
    if [ -f $f ]; then
        basedir=$(basename "$dir")
        echo "Adding: " ${basedir}
        echo ${basedir} >> ${outdir}/test_conditions.txt
    fi
done


# Subset tests for just a few chromosomes in a temp directory
mkdir -p temp/temp_input
for dir in ${inputdir}/*; do 
    f=${dir}/OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv
    if [ -f $f ]; then
        basedir=$(basename "$dir")
        echo $basedir
        mkdir -p temp/temp_input/${basedir}
        cp ${dir}/OPTIM_pcs/base_output__base/cis_nominal1.cis_qtl_pairs.chr1.tsv temp/temp_input/${basedir}
        cp ${dir}/OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv temp/temp_input/${basedir}
    fi
done

