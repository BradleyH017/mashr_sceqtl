# Apply mashr model
# Author: Bradley
# Date: 26/07/24
# module load HGI/softpack/groups/otar2065/seurat5_v2/1
# Applying mashr model to rest of data, per gene

# Load libraries
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
library(rhdf5)
library(zoo)
set.seed(123)

# Define the load function:
h5_mashr_per_gene <- function(f, gene, max.missing=Inf, replace_beta, replace_beta_se) {
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
    message("Dropping ", nrow_start-nrow(betas), " because # NAs > ", 
            max.missing, ". ", nrow(betas), " tests remaining.")

    if(replace_beta  == "median"){
        betas_median <- median(betas, na.rm=TRUE)
        message("Replacing remaining beta NAs with median betas: ", betas_median)
        betas[is.na(betas)] <- betas_median
    } else{
        betas[is.na(betas)] <- as.numeric(replace_beta)
    }

    if(replace_beta_se  == "median"){
        error_median <- median(error, na.rm=TRUE)
        message("Replacing remaining beta NAs with median error: ", error_median)
        error[is.na(error)] <- error_median
    } else {
        error[is.na(error)] <- as.numeric(replace_beta_se)
    }
      
    }
  
    mashData <- mash_set_data(as.matrix(betas), 
                                as.matrix(error))
    return(mashData)
}

# parse args
print("Loading args")
parser <- ArgumentParser()
parser$add_argument("-m", "--model", default = "")
parser$add_argument("-mm", "--max_missing", default = "")
parser$add_argument("-mf", "--merged_file", default = "")
parser$add_argument("-fmb", "--fill_missing_beta", default="./")
parser$add_argument("-fmse", "--fill_missing_se", default="./")
parser$add_argument("-g", "--gene", default="./")
args <- parser$parse_args()
print(args)

# Testing:
# args = args=list(model="results/output/model.rds", max_missing="Inf", merged_file="fastqtl_to_mash_output/merged_test_conditions.h5", fill_missing_beta="0", fill_missing_se="1", gene="ENSG00000000457")

# Load in the data
test_set = h5_mashr_per_gene(args$merged_file, gene = args$gene, max.missing = as.numeric(args$max_missing), replace_beta=args$fill_missing_beta, replace_beta_se=args$fill_missing_se)

# Load in the model
m = readRDS(args$model)

# Apply to this data
m2 = mash(test_set, g=get_fitted_g(m), fixg=TRUE)

# Extract results and save
lfsr <- get_lfsr(m2)
pm <- get_pm(m2)
psd <- get_psd(m2)
groups = colnames(lfsr)
