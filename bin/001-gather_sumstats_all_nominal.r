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
parser$add_argument("-i", "--input_dir", default = "")
parser$add_argument("-cf", "--conditions_file", default = "")
parser$add_argument("-ch", "--chromosomes", default = c(1:22))
parser$add_argument("-s", "--suffix", default = "")
parser$add_argument("-o", "--outdir", default = default="./")
parser$add_argument("-n", "--nrand", default = 10000)
parser$add_argument("-fmb", "--fill_missing_beta", default = 0)
parser$add_argument("-fmse", "--fill_missing_se", default = 10)


# Define combine option
make_mash_dfs <- function(df_list, column, replace_NA) {
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

# Load in the conditions from the input_file
conditions = readLines(args$input_file)

# Load in the results per condition and chromsome
chrs = parser$chrs # NOTE: Adjust this somehow
res = list()
for (c in 1:length(conditions)){
    print(condition)
    per_chr = vector("list", length = length(chrs))
    for(chr in 1:length(chrs)){
        print(chr)
        fname = paste0(inputdir, "/", condition, "/", args$suffix, chr, ".tsv")
        per_chr[[chr]] = read.delim(fname)
    }
    temp = do.call(rbind, per_chr)
    temp$condition = condition
    rownames(temp) = paste(temp$phenotype_id, temp$variant_id, sep = "-")
    res[[condition]] = temp
}

# Then subset the unique tests for value of nrand
if( args$nrand > 0 ){
    all_tests <- c()
    for (df in res) {
        all_tests <- c(all_tests, rownames(df))
    }
    unique_tests <- unique(all_tests)
    # Find max se across all tests
    max_per_cond = lapply(res, function(x) {max(res$slope_se)})
    max_overall <- max(unlist(max_per_cond))
    print(paste0("The maximum SE observed across all tests is: ", max_overall))
    # Subset to nrand tests
    reserve = res
    random_tests <- sample(unique_tests, args$nrand)
    for(r in seq_along(res)){
        res[[r]] = res[[r]][rownames(res[[r]]) %in% random_tests,]
    }
}

# Make the mashr object from this subset
beta_df <- make_mash_dfs(df_list=res, column="slope", replace_NA=args$fill_missing_beta)
if(args$fill_missing_se == "max"){
    se_df <- make_mash_dfs(df_list=res, column="slope_se", replace_NA=max_overall)
} else { 
    se_df <- combine_dataframes(df_list=res, column="slope_se", replace_NA=args$fill_missing_se)
} # From here: https://github.com/stephenslab/mashr/issues/17 somewhat reasonable way to deal with missing data is to set SE to something very large

# Plot the sparsity and missingness
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
ggsave(filename = paste0(args$outdir, "/input/mashr_test_frequency.png"), plot = p, width = 8, height = 6)

# Also plot per condition
library(dplyr)
library(tibble)
count_zeros_per_condition = sapply(beta_df, function(col) sum(col == 0)) %>% 
                                      as.data.frame() %>% 
                                      rownames_to_column(var = "condition") %>%
                                      rename(n_missing_tests = ".")
count_zeros_per_condition$perc_performed_tests = 100*(10000 - count_zeros_per_condition$n_missing_tests) / 10000
med = median(count_zeros_per_condition$perc_performed_tests)
p <- ggplot(count_zeros_per_condition, aes(x = condition, y=perc_performed_tests)) +
  geom_bar(stat="identity") +
  labs(title = "Completeness of 10k test subset across condition",
       x = "Condition",
       y = "Percentage of subset tested (%)") +
  geom_hline(yintercept = med, linetype = "dotted", color = "black") +
  annotate("text", x = Inf, y = med, label = paste0("Median: ", signif(med, 2)), hjust = 1.1, vjust = -0.5, color = "black") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = paste0(args$outdir, "/mashr_test_frequency_per_condition.png"), plot = p, width = 8, height = 6)

# Make the mashr data object
use = list(as.matrix(beta_df), as.matrix(se_df))
names(use) = c("Bhat", "Shat")

# Save the mashr object being used
saveRDS(use, paste0(args$outdir, "/input/sumstats_subset.rds"))



