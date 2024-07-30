# Prep data for mashr
# Author: Bradley
# Date: 26/07/24
# Following: https://stephenslab.github.io/mashr/articles/eQTL_outline.html
# module load HGI/softpack/groups/otar2065/seurat5_v2/1

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
library(rhdf5)
library(zoo)
set.seed(123)


h5_mashr_per_gene <- function(f, gene, max.missing=Inf, replace_beta, replace_beta_se, verbose=F) {
    temp = h5read(f, paste0("/", gene))
    rownames <- temp$rownames
    # add gene
    rownames = paste0(gene, "_", rownames)
    colnames <- temp$colnames
    
    betas <- t(temp$beta)
    rownames(betas) <- rownames
    colnames(betas) <- colnames

    error <- t(temp$se)
    row.names(error) <- rownames
    colnames(error) <- colnames
  
  
  if(any(!complete.cases(betas))){
    nrow_start <- nrow(betas)
    keep <- rowSums(is.na(error)) <= max.missing & rowSums(is.na(betas)) <= max.missing
    error <- error[keep, ]
    betas <- betas[keep, ]
    if(verbose){
      message("Dropping ", nrow_start-nrow(betas), " because # NAs > ", 
              max.missing, ". ", nrow(betas), " tests remaining.")

      message("Have ", signif(100*sum(is.na(betas))/length(betas),2), "% missing values. Filling...")
    }
    if(replace_beta  == "median"){
        betas_median <- median(betas, na.rm=TRUE)
        if(verbose){
          message("Replacing remaining beta NAs with median betas: ", betas_median)
        }
        betas[is.na(betas)] <- betas_median
    } else{
        betas[is.na(betas)] <- as.numeric(replace_beta)
    }

    if(replace_beta_se  == "median"){
        error_median <- median(error, na.rm=TRUE)
        if(verbose){
          message("Replacing remaining beta NAs with median error: ", error_median)
        }
        error[is.na(error)] <- error_median
    } else {
        error[is.na(error)] <- as.numeric(replace_beta_se)
    }
      
    }
  
    mashData <- mash_set_data(as.matrix(betas), 
                                as.matrix(error))
    return(mashData)
}


print("Loading args")
parser <- ArgumentParser()
parser$add_argument("-i", "--input_f", default = "")
parser$add_argument("-mm", "--max_missing", default = "Inf")
parser$add_argument("-rb", "--replace_beta", default = "0")
parser$add_argument("-rb", "--replace_se", default = "1")
parser$add_argument("-qv", "--qval_f", default = "")
parser$add_argument("-ls", "--lfsr_strong", default="0.1")
parser$add_argument("-nr", "--nrand", default="10000")
parser$add_argument("-o", "--outdir", default="./")
parser$add_argument("-m", "--matrices", default="./")
parser$add_argument("-r", "--reference", default="./")
#parser$add_argument("-nc", "--nchunks", default = 10000)
args <- parser$parse_args()
print(args)

# Testing
# args=list(input_h5 = "results/old/chr1_8conds/input/fastqtl_to_mash_output/merged_test_conditions.h5", input_mash = "results/old/chr1_8conds/input/fastqtl_to_mash_output/merged_test_conditions.mash.rds", max_missing="Inf", replace_beta = "0", replace_se = "1", qval_f = "results/input/all_qval.tsv", lfsr_strong="0.1", nrand = "100000", outdir = "results/old/chr1_8conds/output", matrices = "pca,canonical,flash", reference="None")

######### Prepping data ###########
# load in the complete merged dataframe matrix
h5file <- H5Fopen(args$input_h5)
overview = h5ls(h5file)
genes = overview$name[grep("ENSG", overview$name)]
mash_per_gene = lapply(genes, function(x){
    print(x)
    h5_mashr_per_gene(f=args$input_f, g=x, replace_beta="after", replace_beta_se="after")
})

# Combine into a single object
beta_all = lapply(mash_per_gene, function(x){
    return(x$Bhat)
})
beta = do.call(rbind, beta_all)
se_all = lapply(mash_per_gene, function(x){
    return(x$Shat)
})
se = do.call(rbind, se_all)

######## For the prediction of effects across conditions, we want to find the top gene-snp pair per condition ########
# Get all qvalues
all_qval = read.delim(args$qval_f, sep = "\t")
all_qval$gene_variant = paste0(all_qval$phenotype_id, "_", all_qval$variant_id)

# ******* TEMPORRY BUG FIX ******** #
all_qval$gene_variant = gsub("_rs", "_s", all_qval$gene_variant) # UNIQUE TO US, and possibly a bug

# Intersect beta/se, for the variants that were top per condition - in the qval file
beta_pred = beta[rownames(beta) %in% all_qval$gene_variant,]
se_pred = se[rownames(se) %in% all_qval$gene_variant,]
# this is the set of variants x QTLs that we want to predict for

# Use the top hit per gene per condition, or strong from workflow as the 'strong' set
use_strong_from_wf = T
if(use_strong_from_wf){
  obj=readRDS("results/old/chr1_8conds/input/fastqtl_to_mash_output/merged_test_conditions.mash.rds")
} else {
  m.1by1 = mash_1by1(mpred)
  strong.subset = get_significant_results(m.1by1, as.numeric(args$lfsr_strong))
  # Summarise the sparsity here
  genes_strong = unlist(strsplit(rownames(beta_top[strong.subset,]), "\\_"))[c(T,F)]
  genes_count = table(genes_strong)
}

######## Get random subset ########
random.subset = sample(1:nrow(beta), args$nrand)
temp.random = mash_set_data(beta[random.subset,], se[random.subset,])

######### Correlation structure ###########
# Estimating correlation structure from the random tests
Vhat = estimate_null_correlation_simple(temp.random)

######### Set up proper strong/random subsets using the correlation structure ########
data.random = mash_set_data(beta[random.subset,], se[random.subset,], V=Vhat)
if(use_strong_from_wf){
  data.strong = mash_set_data(obj$strong.b, obj$strong.s, V=Vhat)
} else {
  data.strong = mash_set_data(beta[strong.subset,], se[strong.subset,], V=Vhat) 
}

# Update reference if desired
if( args$reference != "None"){
  data.random <- mash_update_data(data.random, ref = args$reference)
  data.strong <- mash_update_data(data.strong, ref = args$reference)
}

######### Derive data-driven covariance matrices ########
want_matrices = unlist(strsplit(args$matrices, ","))
dd_matrices = NULL
if("pca" %in% want_matrices){
  U.pca = cov_pca(data.strong,5)
  dd_matrices = c(U.pca)
}
if("flash" %in% want_matrices){
  U.f = cov_flash(data.strong, factors="nonneg", tag="non_neg")
  if(length(dd_matrices) == 0){
    dd_matrices[[1]] = U.f
  } else {
    dd_matrices = c(dd_matrices, U.f)
  }
}

# Apply extreme deconvolution 
U.ed = cov_ed(data.strong, dd_matrices)

########## Now fit the mash model (data-driven + canonical) to random tests ##########
if("canonical" %in% want_matrices){
  U.c = cov_canonical(data.random)
  Ulist = c(dd_matrices, U.c)
} else {
  Ulist = c(dd_matrices)
}
m = mash(data.random, Ulist = Ulist, outputlevel = 1)

# Save
saveRDS(m, paste0(args$outdir, "/model.rds"))

########## Get posterior summaries #######
if(use_strong_from_wf){
  # Predict for the predicted set - now incorporating correlations
  mpred = mash_set_data(beta_top, se_top, V=Vhat)
  m2 = mash(mpred, g=get_fitted_g(m), fixg=TRUE)
} else {
  # Predict for the strong set
  m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
}

########## Make plots to summarise the model #######
if(file.exists(paste0(args$outdir, "/model_summary")) == F){
  dir.create(paste0(args$outdir, "/model_summary"))
}
plotout = paste0(args$outdir, "/model_summary")

# Plot pairwise sharing:
sharing = get_pairwise_sharing(m2)
melt = reshape2::melt(sharing)
p=ggplot(data = melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_bw() +
  labs(x = "Condition1", y = "Condition2", fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(plotout, "/pairwise_sharing.png"), plot = p, width = 8, height = 6)

# Plot a few very strong hits
for(i in seq(1,5)){
  pdf(file=paste0(plotout, "/top_hit_", i, ".pdf"))
  mash_plot_meta(m2,get_significant_results(m2)[i], thresh=0.05)
  dev.off() 
}

# Save a list of these at varying LFSR thresh
thresh = c(0.01, 0.05, 0.1)
for(t in thresh){
  sig = as.data.frame(get_significant_results(m2, thresh = t))
  sig$effect = rownames(sig)
  rownames(sig) = NULL
  colnames(sig) = c("index", "effect")
  write.table(sig, paste0(plotout, "/predicted_results/sig_effects_any_lfsr_", t, ".txt"), quote=F, row.names=F)
}

# Print/save the log-likelihood
print(get_loglik(m))
write.table(get_loglik(m), paste0(plotout, "/log_likelihood.txt"), col.names=F, quote=F, row.names=F)

# Save the model sharing
write.table(sharing, paste0(plotout, "/pairwise_sharing.txt"))

# plot / save the estimates pi for sharing across covariance matrices
pdf(file=paste0(plotout, "/estimated_pi_covariances.pdf"))
barplot(get_estimated_pi(m),las = 2)
dev.off()

# Make per test results for the strong and random subsets
lfsr <- get_lfsr(m2)
pm <- get_pm(m2)
psd <- get_psd(m2)
groups = colnames(lfsr)
for(i in seq_along(groups)){
  df = data.frame("posterior_means" = pm[,i], "posterior_sd" = psd[,i], "lfsr" = lfsr[,i])
  write.table(df, paste0(plotout, "/predicted_results/", groups[i], ".txt"))
}

# Save a list of uniquely significant effects at varying thresh
for( t in thresh){
  sig_count = apply(lfsr, 1, function(x){sum(x < t)})
  unique = names(sig_count[sig_count == 1])
  write.table(unique, paste0(plotout, "/predicted_results/unique_effects_lfsr_", t, ".txt"), row.names=F, col.names=F, quote=F)
}