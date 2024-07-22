######## Bradley June 2024 #######
# MASHR - COMMON BASELINE AT MEAN. FOLLOWING: https://stephenslab.github.io/mashr/articles/intro_mashbaselinemean.html
# Combining results across tissues x cell-types and looking for effect sharing
# module load HGI/softpack/groups/otar2065/seurat5_v2/1
library('mashr')
library(ggplot2)
library(purrr)
setwd("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/ashg_eqtl")

resdir = "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/ashg_eqtl"
tissues_levels <- list.dirs(resdir, recursive = FALSE, full.names = FALSE)
tissues_levels = tissues_levels[-which(tissues_levels == "temp_out")]
res <- list()


#### TESTING - try for just one tissue / level
tissues_levels = c("ti_label")

# Raw results
nchr=5
for(tl in tissues_levels){ 
    subdirs = list.dirs(paste0(resdir, "/", tl, "/results/TensorQTL_eQTLS"), recursive = FALSE, full.names = T)
    for (dir in subdirs){
        split = unlist(strsplit(dir, "\\/"))
        condition = split[length(split)]
        print(condition)
        per_chr = vector("list", length=nchr)
        if(file.exists(paste0(dir, "/OPTIM_pcs/base_output__base/Cis_eqtls_qval.tsv"))){
          # Only load if was successfully completed
          for(chr in 1:nchr){
            print(chr)
            fname = paste0(dir, "/OPTIM_pcs/base_output__base/cis_nominal1.cis_qtl_pairs.chr", chr, ".tsv")
            print(fname)
            per_chr[[chr]] = read.delim(fname)
          }
          temp = do.call(rbind, per_chr)
        } 
        temp$condition = condition
        rownames(temp) = paste(temp$phenotype_id, temp$variant_id, sep = "-")
        res[[condition]] = temp
      }
    }

# Analyse sparsity in tests
# Number of tested genes, variants and tests (genes x variants) (overall and common)
get_total_common = function(df, col){
  total = length(unique(unlist(lapply(df, function(x) x[,col]))))
  common = length(Reduce(intersect, lapply(df, function(x) x[,col])))
  return(paste0("Total ", col, " tested = ", total, ". Common tested = ", common))
}
genes = get_total_common(res, "phenotype_id")
# "Total phenotype_id tested = 23713. Common tested = 1973"


# Combine together
combine_dataframes <- function(df_list, column, replace_NA) {
  # Initialize an empty dataframe to store the combined data
  combined_df <- data.frame(effect = character(0), stringsAsFactors = FALSE)
  # Loop through the list of dataframes
  for (i in 1:length(df_list)) {
    # Extract the 'effect', 'beta', and 'condition' columns from the current dataframe
    effect_column <- rownames(df_list[[i]])
    beta_column <- df_list[[i]][,column]
    condition <- unique(df_list[[i]]$condition)
    
    # Create a new dataframe for the current condition
    condition_df <- data.frame(effect = effect_column, condition = beta_column, stringsAsFactors = FALSE)
    
    # Rename the 'condition' column to match the condition name
    colnames(condition_df) <- c("effect", condition)
    
    # Merge the current condition's dataframe with the combined dataframe using 'effect' as the key
    combined_df <- merge(combined_df, condition_df, by = "effect", all = TRUE)
  }
  # Set 'effect' column as row names
  rownames(combined_df) <- combined_df$effect
  combined_df$effect <- NULL
  # Replace NA values with 0
  combined_df[is.na(combined_df)] <- replace_NA
  return(combined_df)
}

### Subset to 10k random tests from total sets
# Get unique vector of all tests
all_tests <- c()
for (df in res) {
  all_tests <- c(all_tests, rownames(df))
}
unique_tests <- unique(all_tests)
# Subset to top 10k tests from first condition
reserve = res
set.seed(123)
random_tests <- sample(unique_tests, 10000)
for(r in seq_along(res)){
  res[[r]] = res[[r]][rownames(res[[r]]) %in% random_tests,]
}

# Call to generate the beta and se matrices 
beta_df <- combine_dataframes(df_list=res, column="slope", replace_NA=0)
# Find the percentage missing:
count_zeros <- sum(sapply(beta_df, function(col) sum(col == 0)))
perc_missing = count_zeros/(ncol(beta_df)*nrow(beta_df))
se_df <- combine_dataframes(df_list=res, column="slope_se", replace_NA=max(unlist(lapply(res, function(x){max(x$slope_se)})))) # From here: https://github.com/stephenslab/mashr/issues/17 somewhat reasonable way to deal with missing data is to set SE to something very large

# Have a look at sparsity of tests
# Draw the distribution of test performance
non_zero_counts <- rowSums(beta_df != 0)
df_non_zero_counts <- data.frame(non_zero_counts = non_zero_counts)
p <- ggplot(df_non_zero_counts, aes(x = non_zero_counts)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = paste0("Distribution of Non-Zero Counts in test subset \nTotal missingness = ", signif(perc_missing, 2), "%"),
       x = "Number of Non-Zero Values",
       y = "Frequency") +
  theme_bw()

# Save the plot to a file using ggsave
ggsave(filename = "temp_out/mashr_test_frequency.png", plot = p, width = 8, height = 6)

# Also plot per condition
library(dplyr)
library(tibble)
count_zeros_per_condition = sapply(beta_df, function(col) sum(col == 0)) %>% 
                                      as.data.frame() %>% 
                                      rownames_to_column(var = "condition") %>%
                                      rename(n_missing_tests = ".")
count_zeros_per_condition$perc_performed_tests = 100*(10000 - count_zeros_per_condition$n_missing_tests) / 10000
p <- ggplot(count_zeros_per_condition, aes(x = condition, y=perc_performed_tests)) +
  geom_bar(stat="identity") +
  labs(title = "Percentage of subset tested in each condition",
       x = "Condition",
       y = "Performed test (%)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = "temp_out/mashr_test_frequency_per_condition.png", plot = p, width = 8, height = 6)

# Make the mashr data object
use = list(as.matrix(beta_df), as.matrix(se_df))
names(use) = c("Bhat", "Shat")

# Randomly sample for 10000 tests
data = mash_set_data(use$Bhat, use$Shat)

# Update to have the reference as the mean
data.L = mash_update_data(data, ref = 'mean')

# Derive the canonical matrices
U.c = cov_canonical(data.L)

# Derive the data-driven covariances
m.1by1 = mash_1by1(data.L)
strong = get_significant_results(m.1by1)
U.pca = cov_pca(data.L, npc=2, subset=strong)
U.ed = cov_ed(data.L, U.pca, subset=strong)

# Now fit the mashr model
m = mash(data.L, c(U.c, U.ed, U.pca), algorithm.version = 'R')

# Apply the mashr model to the total input
# Data must be in the form of original input
res_test = lapply(reserve, function(x){tail(x, 1000)})
beta_test <- combine_dataframes(df_list=res_test, column="slope", replace_NA=0)
se_test <- combine_dataframes(df_list=res_test, column="slope_se", replace_NA=max(unlist(lapply(res, function(x){max(x$slope_se)})))) # From here: https://github.com/stephenslab/mashr/issues/17 somewhat reasonable way to deal with missing data is to set SE to something very large
use_test = list(as.matrix(beta_test), as.matrix(se_test))
names(use_test) = c("Bhat", "Shat")
data_test = mash_set_data(use_test$Bhat, use_test$Shat)
r1 = mash(data_test, g=get_fitted_g(m), fixg=TRUE)