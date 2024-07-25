# Prep data for mashr
# Author: Bradley
# Date: 22/07/24
# Following: https://stephenslab.github.io/mashr/articles/eQTL_outline.html
# module load HGI/softpack/groups/otar2065/seurat5_v2/1

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
set.seed(123)

print("Loading args")
parser <- ArgumentParser()
parser$add_argument("-i", "--input_f", default = "")
parser$add_argument("-o", "--outdir", default="./")
parser$add_argument("-m", "--matrices", default="./")
parser$add_argument("-r", "--reference", default="./")
#parser$add_argument("-nc", "--nchunks", default = 10000)
args <- parser$parse_args()
print(args)

# Testing
# args=list(input_f = "results/input_chr1_4/fastqtl_to_mash_output/merged_test_conditions.mash.rds", outdir = "results/output", matrices = "pca,canonical,flash", reference=NULL)

# Load in the mashr object
obj <- readRDS(args$input_f)

######### Correlation structure ###########
# Estimating correlation structure from the random tests
Vhat = estimate_null_correlation_simple(mash_set_data(obj$random.b, obj$random.s))

######### Set up proper stron/random subsets using the correlation structure ########
data.random = mash_set_data(obj$random.b, obj$random.s, V=Vhat)
data.strong = mash_set_data(obj$strong.b, obj$strong.s, V=Vhat)

# Update reference if desired
if(!is.null(args$reference)){
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
saveRDS(m, "results/output/model.rds")

########## Get posterior summaries #######
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

########## Make plots to summarise the model #######
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
  mash_plot_meta(m2,get_significant_results(m2)[i])
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