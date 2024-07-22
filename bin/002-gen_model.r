# Apply fit mashr model to data
# Author: Bradley
# Date: 22/07/24

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(argparse)
set.seed(123)

message("Loading arguments...")
parser <- ArgumentParser()
parser$add_argument("-i", "--input_file", default = "")
parser$add_argument("-o", "--outdir", default = default="./")
parser$add_argument("-r", "--reference", default = default=NULL)
parser$add_argument("-l", "--lfsr", default = 0.1)
args <- parser$parse_args()

# Load the input data
use = readRDS(args$input_file)

# Make the mash object
data = mash_set_data(use$Bhat, use$Shat)

# Define reference (if desired)
if(!is.null(variable)){
    data = mash_update_data(data, ref = args$reference)
}

# Derive the canonical matrices
U.c = cov_canonical(data)

# Derive the data-driven covariances
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1, thresh = args$lfsr)
U.pca = cov_pca(data, npc=5, subset=strong) 
U.ed = cov_ed(data, U.pca, subset=strong)

# Now fit the mashr model
m = mash(data.L, c(U.c, U.ed, U.pca), algorithm.version = 'R')

# Save this
saveRDS(m, paste0(outdir, "/mash_model/mash.rds"))

# Get pairwise sharing
sharing = get_pairwise_sharing(m)
melt = reshape2::melt(sharing)
p=ggplot(data = melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_bw() +
  labs(x = "Condition1", y = "Condition2", fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(outdir, "/mash_model/pairwise_sharing.png"), plot = p, width = 14, height = 12)

# Plot 5 really strong hit
for(c in c(1:5)){
    pdf(paste0(outdir, "/mash_model/strong_hit", c, ".png"))
    mash_plot_meta(m,get_significant_results(m)[c])
    dev.off()
}

# Number of significant results:
nsig = (get_significant_results(m, thresh = 0.05)) # default lfsr = 0.05
print(paste0("The number of significant hits at LFSR < 0.05 = ", length(nsig)))

# Summarise model fit
print("The overall model fit is")
print(get_loglik(m))


