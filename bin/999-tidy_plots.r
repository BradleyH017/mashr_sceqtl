# Load packages
library('mashr')
library(ggplot2)
library(purrr)
library(flashier)
library(argparse)
library(rhdf5)
library(zoo)
library(ggdendro)
library(dendextend)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(dplyr)
library(tidyr)
set.seed(123)

####### 1.  Plot sharing
sharing = read.delim("results/output/model_summary/pairwise_sharing.txt", sep = " ")
rownames(sharing) = gsub("merged_", "", rownames(sharing))
rownames(sharing) = gsub("-dMean.tsv", "", rownames(sharing))
colnames(sharing) = gsub("merged_", "", colnames(sharing))
colnames(sharing) = gsub(".dMean.tsv", "", colnames(sharing))
colnames(sharing) = gsub("\\.", "\\-", colnames(sharing))

# Plot with clustered rows and columns
row_hc <- hclust(dist(sharing, method = "euclidean"), method = "complete")
col_hc <- hclust(dist(t(sharing), method = "euclidean"), method = "complete")

ppi=400
png("results/output/model_summary/tidy_sharing_matrix.png", res=ppi, width=8*ppi, height=8*ppi)
Heatmap(as.matrix(sharing),
        cluster_rows = row_hc,
        cluster_columns = col_hc,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        col = colorRamp2(c(0,1), c("white", "red")),
        name = "Sharing",
        column_title = "Annotation 1",
        row_title = "Annotation 2",
        heatmap_legend_param = list(title = "Sharing"),
        #row_dend_side = "right", column_dend_side = "bottom",
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),  
        row_dend_width = unit(8, "mm"), 
        column_dend_height = unit(8, "mm"))

dev.off()

# Find the median sharing within category vs outside
conv = read.csv("/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/sc-eqtl_meta/Anderson_eqtl_generation/input/new_TI_v7_atlas_to_Franke.csv")
to_add = conv[,c("manual_annotation_retired_machine", "category")]
sharing$manual_annotation_retired_machine = unlist(strsplit(rownames(sharing), "\\-"))[c(F,T)]
sharing$query_annot = rownames(sharing)
melt = melt(sharing)
colnames(melt)[which(colnames(melt) == "variable")] = "ref_annot"
melt = merge(melt, to_add, by="manual_annotation_retired_machine")
colnames(melt)[which(colnames(melt) == "category")] = "query_cat"

# Plot within vs outside medians
cats = unique(melt$query_cat)
inout = vector("list", length = length(cats))
for(c in cats){
    print(c)
    inside_labs = unique(melt[melt$query_cat == c,]$query_annot)
    outside_labs = setdiff(melt$ref_annot, inside_labs)
    inside_sharing = median(melt[melt$query_annot %in% inside_labs & melt$ref_annot %in% inside_labs & melt$query_annot != melt$ref_annot, ]$value)
    outside_sharing = median(melt[melt$query_annot %in% inside_labs & melt$ref_annot %in% outside_labs & melt$query_annot != melt$ref_annot, ]$value)
    temp = data.frame("category" = c, "median_in_category_sharing" = inside_sharing, "median_outside_category_sharing" = outside_sharing)
    inout[[which(cats == c)]] = temp
}
inout_all = do.call(rbind, inout)

# Plot this (paired, lines between the two)
inout_long <- inout_all %>%
  pivot_longer(cols = c(median_in_category_sharing, median_outside_category_sharing),
               names_to = "sharing_type",
               values_to = "value")

inout_long$sharing_type <- factor(inout_long$sharing_type, 
                                  levels = c("median_in_category_sharing", "median_outside_category_sharing"),
                                  labels = c("In Category", "Outside Category"))


p = ggplot(inout_long, aes(x = sharing_type, y = value)) +
    geom_boxplot() + 
    geom_point(aes(color = category), size = 3) +  
    geom_line(aes(group = category)) + 
    theme_bw() +
    labs(x = "Sharing Type", y = "Median sharing value", title = "Comparison of sharing vs transcriptional similarity") +
    theme(legend.title = element_blank())  
ggsave(filename = "results/output/model_summary/Sharing_comparison_in_vs_out.png", plot = p, width = 8, height = 6)

# Look at b cell plasma
bcp = unique(melt[melt$query_cat == "Plasma B cell",]$query_annot)
sharing[rownames(sharing) %in% bcp, colnames(sharing) %in% bcp]


# Check a specific result
resdir="results/output/model_summary/predicted_results"
resf = list.files(resdir)
resf = resf[grep("machine", resf)]
gene = "ENSG00000108175" # ZMIZ1
system(sprintf("mkdir -p %s/temp && for f in %s/*machine*; do echo $f; grep '%s' $f > %s/temp/$(basename $f)_%s.txt; done", 
               resdir, resdir, gene, resdir, gene))

filtresdir = paste0(resdir, "/temp")
filtresdirf = list.files(filtresdir)
filtresdirf = filtresdirf[grep(gene, filtresdirf)]
res_all = vector("list", length = length(filtresdirf))
for(f in filtresdirf){
    print(f)
    temp = read.delim(paste0(filtresdir, "/", f), header=F, sep = " ")
    temp$f = f
    res_all[[which(filtresdirf == f)]] = temp
}
resdf = do.call(rbind, res_all)
colnames(resdf) = c("effect", "posterior_means", "posterior_sd", "lfsr", "source")
resdf$annotation = gsub(paste0("-dMean_mashr.txt_", gene, ".txt"), "", resdf$source)


# Custom forest plot for specific gene and variants
var = "s1250556"
gene = "ENSG00000108175"
resdf_sub = resdf[grep(gene, resdf$effect), ]
resdf_sub = resdf_sub[grep(var, resdf_sub$effect),]

resdf_sub$significance <- ifelse(resdf_sub$lfsr < 0.0001, "****",
                         ifelse(resdf_sub$lfsr < 0.001, "***",
                         ifelse(resdf_sub$lfsr < 0.01, "**",
                         ifelse(resdf_sub$lfsr < 0.05, "*", ""))))


ggplot(resdf_sub, aes(x = posterior_means, y = annotation)) +
  geom_point() +  
  geom_errorbarh(aes(xmin = posterior_means - posterior_sd, 
                     xmax = posterior_means + posterior_sd), height = 0.2) + 
  geom_text(aes(label = significance), vjust = -0.15, hjust = -0.5, size = 2.5) + 
  theme_bw() +
  geom_vline(xintercept=0, lty="dashed", col="grey") + 
  labs(x = "Posterior Means", y = "Annotation", title = unique(resdf_sub$effect)) +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 10))
dev.off()