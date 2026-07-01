#create a violin plot with vst normalized data and compare the counts among genes of interest

create_violin_plot_compare_genes_vst <- function(dds, genes, intgroup,return_data=F, select_group = NULL, group_labs = NULL, ncol = 4) {
  # Load necessary libraries
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(tidyr)
  
  # Check if intgroup is in colData(dds)
  if (!(intgroup %in% colnames(colData(dds)))) {
    stop(paste("The intgroup variable", intgroup, "is not in colData(dds)."))
  }
  
  # Get the VST-transformed counts
  vst_data <- assay(vst(dds))[genes, , drop=FALSE]
  vst_data <- as.data.frame(vst_data)
  vst_data$genes <- rownames(vst_data)
  
  # Reshape data to long format
  vst_data_long <- vst_data %>%
    pivot_longer(cols = -genes, names_to = "sample", values_to = "vst_count")
  
  # Get intgroup variable from colData
  coldata <- as.data.frame(colData(dds))
  coldata$sample <- rownames(coldata)
  
  # Merge vst_data_long with coldata
  vst_data_long <- merge(vst_data_long, coldata, by="sample")
  
  # Ensure vst_data_long contains intgroup
  if (!(intgroup %in% colnames(vst_data_long))) {
    stop(paste("The intgroup variable", intgroup, "is not in vst_data_long."))
  }
  
  # Subset vst_data_long if select_group is specified
  if (!is.null(select_group)) {
    vst_data_long <- vst_data_long[vst_data_long[[intgroup]] %in% select_group, ]
    vst_data_long[[intgroup]] <- droplevels(factor(vst_data_long[[intgroup]]))
  }
  
  # If group labels are provided, rename the factor levels
  if (!is.null(group_labs)) {
    levels(vst_data_long[[intgroup]]) <- group_labs
  }
  
  # Convert genes and intgroup to factors
  vst_data_long$genes <- as.factor(vst_data_long$genes)
  vst_data_long$genes <- factor(vst_data_long$genes, 
                               levels = genes)
  vst_data_long[[intgroup]] <- as.factor(vst_data_long[[intgroup]])

  
  # Perform pairwise comparisons between genes within each group
  
    if(length(genes) == 1){
      stats_res <- vst_data_long %>%
        group_by(genes) %>%
        pairwise_t_test(
          vst_count ~ KNN_cluster,
          p.adjust.method = "BH"
        )%>%
        add_significance("p.adj") %>%
        add_xy_position(x = "KNN_cluster", step.increase = 1)
    } else {
      stats_res <- vst_data_long %>%
      dplyr::group_by(!!sym(intgroup)) %>%
        pairwise_t_test(
          vst_count ~ genes,
          p.adjust.method = "BH"
        ) %>%
        add_significance("p.adj") %>%
        add_xy_position(x = "genes")
      }
  
  # Create a violin plot comparing genes
  if(length(genes) > 1){
    p <- ggplot(vst_data_long, aes(x = genes, y = vst_count, fill = genes)) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.2) +
      facet_wrap(as.formula(paste("~", intgroup)), ncol = ncol) +
      labs(title = "Violin Plot Comparing genes", y = "VST-Transformed Counts") +
      theme_classic() +
      theme(legend.position = "bottom")+ stat_pvalue_manual(
        stats_res,
        label = "p.adj.signif",
        inherit.aes = FALSE,size = 7)
  } else {
    p <- ggplot(vst_data_long, aes(x = KNN_cluster, y = vst_count, fill = genes, group=KNN_cluster)) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.2) +
      labs(title = "Violin Plot Comparing genes", y = paste0("VST-Transformed Counts of ", genes)) +
      theme_classic() +
      theme(legend.position = "bottom")+ stat_pvalue_manual(
        stats_res,
        label = "p.adj",
        inherit.aes = FALSE,size = 7)
  }
  

  if(return_data){
    return(list(
      plot = p, 
      data = vst_data_long))
  } else {
    return(p)
  }
}
