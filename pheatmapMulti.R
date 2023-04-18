
#A function to plot multiple categories

pheatmapMulti <- function(deseq2_outputs_list, normcounts, metafile){
  # Get_normalized counts
  normcounts_df <- data.frame(normcounts)
  normcounts_df["genes"] <- rownames(normcounts_df)
  # Get metafile
  metafile <- read.table("metafile_corrected.txt")
  colnames(metafile) <- c("genes", "features")
  featured_genes <- metafile$genes
  normcounts_df_meta <- merge(normcounts_df, metafile, by="genes")
  # Prepare master gene list for all DE genes
  master_deseq2_outputs_genes <- unique(unlist(lapply(deseq2_outputs_list, function(x) rownames(x))))
  # Get list of all DE genes corresponding to Samples
  subList_deseq2_outputs <- lapply(deseq2_outputs_list, function(x) (ifelse(master_deseq2_outputs_genes%in%rownames(x) ==TRUE, yes = "DE", no = "nonDE")))
  # Combined Samples
  subDF_deseq2_outputs <- do.call(cbind, subList_deseq2_outputs)
  rownames(subDF_deseq2_outputs) <- master_deseq2_outputs_genes
  colnames(subDF_deseq2_outputs) <- gsub("_results_0.05","",gsub("deseq2_pairwise_", "", colnames(subDF_deseq2_outputs)))
  subDF_deseq2_outputs <- data.frame(subDF_deseq2_outputs)
  subDF_deseq2_outputs["genes"] <- rownames(subDF_deseq2_outputs)
  # Extract DE genes
  normcounts_df_meta_de <- merge(normcounts_df_meta, subDF_deseq2_outputs, by="genes")
  rownames(normcounts_df_meta_de) <- normcounts_df_meta_de$genes
  # Prepare pheatmap annotations
  my_gene_col <- normcounts_df_meta_de[,c((length(normcounts_df_meta_de)-(length(subDF_deseq2_outputs)-1)):(length(normcounts_df_meta_de)))]
  my_colour = list()
  for (categories in colnames(my_gene_col)){
    print(categories)
    my_colour[[categories]] <- c(DE = "yellow", nonDE = "blue")
  }
  # "features" column need to be removed as it is not about DE
  my_colour$features <- NULL
  # Feature free df: Just cleanup
  normcounts_df_meta_de_featured <- normcounts_df_meta_de[,c(2:(length(normcounts_df_meta_de)-length(subDF_deseq2_outputs)))]
  breaksListn <- seq(-4, 4, by = 0.01)
  # Plot heatmap
  pheatmap(normcounts_df_meta_de_featured,color = colorRampPalette(c("red", "black", "green"))(length(breaksListn)), breaks = breaksListn,fontsize = 6,clustering_distance_cols = "euclidean",cluster_rows = T,cluster_cols = T,clustering_method = "ward.D", border_color = "black",show_rownames = T, show_colnames = T, annotation_row = my_gene_col,annotation_colors = my_colour)
}


#Run pheamtpamulti function as follows:
#prepare list of deseq2 objects
deseq2_outputs_list <- list(deseq2_pairwise_Mutant_CYC_TP1_Wildtype_CYC_WP1_results_0.05= deseq2_pairwise_Mutant_CYC_TP1_Wildtype_CYC_WP1_results_0.05, deseq2_pairwise_Mutant_CYC_TP2_Wildtype_CYC_WP2_results_0.05=deseq2_pairwise_Mutant_CYC_TP2_Wildtype_CYC_WP2_results_0.05, deseq2_pairwise_Mutant_CYC_TP3_Wildtype_CYC_WP3_results_0.05=deseq2_pairwise_Mutant_CYC_TP3_Wildtype_CYC_WP3_results_0.05)
#lets say you want z-score
z_Tason_ddsONormcounts= scale(t(ason_ddsONormcounts), center = TRUE, scale = TRUE)
z_ason_ddsONormcounts <- data.frame(t(z_Tason_ddsONormcounts))
head(z_ason_ddsONormcounts)

pheatmapMulti(deseq2_outputs_list, z_ason_ddsONormcounts, metafile)
