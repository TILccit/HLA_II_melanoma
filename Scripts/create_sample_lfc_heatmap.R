#function to create sample log fold change heatmap to visualize 
#transcriptomic changes among the clusters

create_sample_lfc_heatmap <- function(se_object, genes, high_expression_ids, low_expression_ids) {
  # Extract the assay data
  data <- as.data.frame(assay(se_object))
  
  # Ensure only the selected genes are used and handle missing genes
  missing_genes <- vector()
  for (g in genes) {
    if (!g %in% rownames(data)) {
      print(paste0("removing ", g, " from genes as it is not present in the expression data"))
      missing_genes <- c(missing_genes, g)
    }
  }
  genes <- genes[!genes %in% missing_genes]
  data_subset <- data[genes, ]
  
  # Create sample labels based on high and low expression IDs
  sample_labels <- ifelse(colnames(data_subset) %in% high_expression_ids, "HLA_high",
                          ifelse(colnames(data_subset) %in% low_expression_ids, "HLA_low", "HLA_mid"))
  sample_labels <- factor(sample_labels, levels = c("HLA_high", "HLA_low", "HLA_mid"))
  
  # Calculate log fold change between each sample
  # logFC = log2(high_sample / low_sample)
  logfc_matrix <- matrix(nrow = length(genes), ncol = ncol(data_subset))
  rownames(logfc_matrix) <- genes
  colnames(logfc_matrix) <- colnames(data_subset)
  
  for (i in 1:ncol(data_subset)) {
    # Log fold change for each sample against the mean of the opposite group
    if (colnames(data_subset)[i] %in% high_expression_ids) {
      mean_low <- rowMeans(data_subset[, low_expression_ids])
      logfc_matrix[, i] <- log2(data_subset[, i] / mean_low)
    } else if (colnames(data_subset)[i] %in% low_expression_ids) {
      mean_high <- rowMeans(data_subset[, high_expression_ids])
      logfc_matrix[, i] <- log2(data_subset[, i] / mean_high)
    } else {
      logfc_matrix[, i] <- NA  # Handle "mid" or other cases as NA
    }
  }
  
  # Create an annotation data frame
  annotation_col <- data.frame(Group = sample_labels)
  rownames(annotation_col) <- colnames(data_subset)
  
  # Define colors for the annotation
  annotation_colors <- list(Group = c(HLA_high = "red", HLA_low = "blue"))
  
  library(pheatmap)
  
  # Generate the heatmap for log fold change in each sample
  ph <- pheatmap(t(logfc_matrix),
           scale = "none",  # Log fold change doesn't need further scaling
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           annotation_row = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = F,
           show_colnames = F,
           main = "LFC Heatmap of DE Genes")
  legend_index <- which(ph$gtable$layout$name == "legend")
  # add legend title
  library(grid)
  grid.text("Log2 fold change", x = unit(0.85, "npc"), y = unit(0.87, "npc"), just = "left", gp = gpar(fontsize = 10, fontface="bold"))
return(ph)
}
