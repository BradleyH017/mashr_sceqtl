# Prep data for mashr
# Author: Bradley
# Date: 22/07/24
# Following: https://stephenslab.github.io/mashr/articles/eQTL_outline.html

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(argparse)
set.seed(123)

message("Loading arguments...")
parser <- ArgumentParser()
parser$add_argument("-i", "--input_f", default = "")
parser$add_argument("-o", "--outdir", default="./")
parser$add_argument("-nc", "--nchunks", default = 10000)
args <- parser$parse_args()
print(args)

# Testing
# args=list(input_f = "results/input/fastqtl_to_mash_output/merged_test_conditions.mash.rds", outdir = "results/output", nchunks = 1000)


# Load in the mashr object
obj <- readRDS(args$input_f)

######### Strong ###########
# First load in the strong results for each condition
print("****** Generating data-driven covariance matrices from the strong results ******")
m.1by1.strong = mash_1by1(mash_set_data(obj$strong.b, obj$strong.s))
strong.subset = get_significant_results(m.1by1.strong,0.05)

######### Correlation structure ###########
# Now estimating correlation structure from the random tests
Vhat = estimate_null_correlation_simple(mash_set_data(obj$random.b, obj$random.s))

######### Set up proper stron/random subsets using the correlation structure ########
data.random = mash_set_data(obj$random.b, obj$random.s, V=Vhat)
data.strong = mash_set_data(obj$strong.b[strong.subset,], obj$strong.s[strong.subset,], V=Vhat)

######### Derive data-driven covariance matrices ########
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)

########## Now fit the mash model (data-driven + canonical) to random tests ##########
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

########## Get posterior summaries #######
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

########## Make plots to summarise the model #######
plotout = paste0(args$outdir, "/model")

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
pdf(file=paste0(plotout, "/top_hit.pdf"))
mash_plot_meta(m2,get_significant_results(m2)[1])
dev.off()

