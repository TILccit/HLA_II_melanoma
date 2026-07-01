#function to perform cluster analysis in RNASeq data

perform_cluster_analysis <- function(expression_matrix, genes, annotation = NULL, perplexity = ceiling(0.1149425 * num_samples), clusters = NULL) {
  # Extract sample annotations if expression_matrix is DESeqTransform
  if (class(expression_matrix) == "DESeqTransform") {
    sample_annotations <- as.data.frame(colData(expression_matrix))
    expression_matrix <- as.data.frame(assay(expression_matrix))
  } else {
    stop("expression_matrix must be a DESeqTransform object.")
  }
  
  # Validate and extract the annotation
  if (!is.null(annotation)) {
    if (length(annotation) != 1) {
      stop("Only one annotation is allowed.")
    }
    # Check that the annotation is in sample_annotations
    if (!annotation %in% colnames(sample_annotations)) {
      stop(paste("The annotation is not present in the sample annotations:", annotation))
    }
    # Extract the annotation
    selected_annotation <- sample_annotations[[annotation]]
    # Ensure that the annotation is a factor
    selected_annotation <- as.factor(as.character(selected_annotation))
    names(selected_annotation) <- rownames(sample_annotations)
  }
  
  # Filter the expression matrix to include only specified genes
  missing_genes <- vector()
  for (g in genes) {
    if (!g %in% rownames(expression_matrix)) {
      print(paste0("Removing ", g, " from genes as it is not present in the expression data"))
      missing_genes <- c(missing_genes, g)
    }
  }
  genes <- genes[!genes %in% missing_genes]
  
  filtered_data <- expression_matrix[genes, ]
  filtered_data <- as.data.frame(t(filtered_data))
  
  # Calculating distances between samples
  distances <- stats::dist(filtered_data, method = "euclidean")
  
  # Performing hierarchical clustering
  hc <- hclust(distances, method = "ward.D2")
  
  # Cutting the dendrogram to create 4 clusters
  clusters_HC <- as.factor(cutree(hc, k = clusters))
  data_for_cluster_labels <- cbind(filtered_data, cluster_id = clusters_HC)
  cluster_means <- data_for_cluster_labels %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(mean_ciita = mean(CIITA, na.rm = TRUE))
  if (clusters == 2){
    cluster_order <- cluster_means %>%
      arrange(mean_ciita) %>%
      mutate(label = c("HLA_low", "HLA_high"))
  } else if (clusters == 4){
    cluster_order <- cluster_means %>%
      arrange(mean_ciita) %>%
      mutate(label = c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))
  } else if (clusters != 2 | 4){
    print("Please, select 2 or 4 clusters")
  }
  
  cluster_label_mapping <- setNames(cluster_order$label, cluster_order$cluster_id)
  levels(clusters_HC) <- cluster_label_mapping[levels(clusters_HC)]
  
  # Run t-SNE
  num_samples <- nrow(filtered_data)
  tsne_results <- Rtsne(filtered_data, dims = 2, perplexity = perplexity, verbose = TRUE)
  
  # The output coordinates and clustering with k-means
  tsne_data <- as.data.frame(tsne_results$Y)
  rownames(tsne_data) <- rownames(filtered_data)
  set.seed(1234)
  kmeans_result <- kmeans(tsne_data, centers = clusters)  # 4 centers = 4 clusters
  data_for_cluster_labels <- cbind(filtered_data, cluster = kmeans_result$cluster)
  cluster_means <- data_for_cluster_labels %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(mean_ciita = mean(CIITA, na.rm = TRUE))
  if (clusters == 2){
    cluster_order <- cluster_means %>%
      arrange(mean_ciita) %>%
      mutate(label = c("HLA_low", "HLA_high"))
  } else if (clusters == 4){
    cluster_order <- cluster_means %>%
      arrange(mean_ciita) %>%
      mutate(label = c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))
  } else if (clusters != 2 | 4){
    print("Please, select 2 or 4 clusters")
  }
  cluster_label_mapping <- setNames(cluster_order$label, cluster_order$cluster)
  tsne_data$cluster <- factor(kmeans_result$cluster)
  levels(tsne_data$cluster) <- cluster_label_mapping[levels(tsne_data$cluster)]
  
  # Prepare cluster labels for heatmap
  cluster_labels <- data.frame(KNN_cluster = tsne_data$cluster, HC_cluster = clusters_HC)
  rownames(cluster_labels) <- rownames(filtered_data)
  
  ### ADD HERE ADDITIONAL OPTIONAL ANNOTATION FOR TSNE ###
  # Add annotation to tsne_data
  if (!is.null(annotation)) {
    tsne_data[[annotation]] <- selected_annotation[rownames(tsne_data)]
  }
  
  # Define colors for clusters
  if (clusters == 2){
    annotation_colors <- list(
      KNN_cluster = setNames(c("#0384fc", "red"), c("HLA_low", "HLA_high")),
      HC_cluster = setNames(c("#0384fc", "red"), c("HLA_low", "HLA_high"))
    )
  } else if (clusters == 4){
    annotation_colors <- list(
      KNN_cluster = setNames(c("#0384fc", "darkblue", "#8b0000", "red"), c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high")),
      HC_cluster = setNames(c("#0384fc", "darkblue", "#8b0000", "red"), c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))
    )
  }
  
  
  ### ADD HERE ADDITIONAL OPTIONAL ANNOTATION FOR PHEATMAP ###
  ### Annotation to be added to the cluster_labels dataframe ###
  # Add annotation to cluster_labels
  if (!is.null(annotation)) {
    cluster_labels[[annotation]] <- selected_annotation[rownames(cluster_labels)]
    # Generate colors for the new annotation
    levels_ann <- levels(cluster_labels[[annotation]])
    num_levels <- length(levels_ann)
    colors_ann <- setNames(colorRampPalette(brewer.pal(8, "Spectral"))(num_levels), levels_ann)
    annotation_colors[[annotation]] <- colors_ann
  }
  
  # Plotting t-SNE
  tSNE_plot <- ggplot(tsne_data, aes(x = V1, y = V2, color = cluster)) +
    geom_point(alpha = 0.5, size=3) +
    labs(title = "t-SNE and k-Means Clustering Results", x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal()+ scale_color_manual(values = annotation_colors[["KNN_cluster"]])
  
  if (!is.null(annotation)) {
    tSNE_plot <- tSNE_plot + aes_string(shape = annotation)
  }
  
  # Generating the heatmap with annotations
  heatmap <- pheatmap(filtered_data,
                      annotation_row = cluster_labels,
                      annotation_colors = annotation_colors,
                      color = colorRampPalette(c("blue", "white", "red"))(100),
                      scale = "column", 
                      show_rownames = F)
  
  high_expression_ids <- rownames(subset(cluster_labels, KNN_cluster == "HLA_high"))
  low_expression_ids <- rownames(subset(cluster_labels, KNN_cluster == "HLA_low"))
  if (clusters > 2){
    mid_high_expression_ids <- rownames(subset(cluster_labels, KNN_cluster == "HLA_mid_high"))
    mid_low_expression_ids <- rownames(subset(cluster_labels, KNN_cluster == "HLA_mid_low"))
  }
  
  return(list(
    tSNE_plot = tSNE_plot,
    heatmap = heatmap,
    cluster_labels = cluster_labels,
    high_expression_ids = high_expression_ids,
    mid_high_expression_ids = mid_high_expression_ids,
    mid_low_expression_ids = mid_low_expression_ids,
    low_expression_ids = low_expression_ids
  ))
}