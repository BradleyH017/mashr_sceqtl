# Generating input conditions from TensorQTL output
inputdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/ashg_eqtl/rectum_label/results/TensorQTL_eQTLS
outdir=/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/mashr_sceqtl/temp/temp_out

rm ${outdir}/test_conditions.txt
for dir in ${inputdir}/*; do
    f=${dir}/OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv
    if [ -f $f ]; then
        basedir=$(basename "$dir")
        echo "Adding: " ${basedir}
        echo ${basedir} >> ${outdir}/test_conditions.txt
    fi
done