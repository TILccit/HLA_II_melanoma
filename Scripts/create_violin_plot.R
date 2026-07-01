#create violin plot function
library(rstatix)
library(ggpubr)

create_violin_plot <- function(dds, gene, intgroup, ref_group = "HLA_low", select_group = c("HLA_low","HLA_high"), group_labs = NULL) {
  # Check if intgroup is in colData(dds)
  if (!(intgroup %in% colnames(colData(dds)))) {
    stop(paste("The intgroup variable", intgroup, "is not in colData(dds)."))
  }
  
  # Get the normalized counts for the gene of interest
  counts_data <- plotCounts(dds, gene = gene, intgroup = intgroup, returnData = TRUE)
  
  # Ensure counts_data contains intgroup
  if (!(intgroup %in% colnames(counts_data))) {
    stop(paste("The intgroup variable", intgroup, "is not in counts_data."))
  }
  
  # Print counts_data columns for debugging
  print("Columns in counts_data:")
  print(colnames(counts_data))
  
  #to be used only when comparing high vs low
  if (!(is.null(select_group))){
    counts_data <- subset(counts_data, KNN_cluster %in% select_group)
    counts_data$KNN_cluster <- droplevels(counts_data$KNN_cluster)
    #extract p value from deseq results
    diff_expr <- as.data.frame(results(dds, cooksCutoff = T,name = paste0("KNN_cluster_",select_group[2],"_vs_",select_group[1])))
    pval_gene <- diff_expr[gene,]
    pval_gene$group1 <- select_group[1]
    pval_gene$group2 <- select_group[2]
    # if (!(is.null(group_labs))){
    #   counts_data
    #TO BE ADDED A PART OF THE FUNCTION TO CHANGE THE NAME OF THE LEVEL $INTGROUP ACCORDING TO GROUP NAMES
    pval_gene$xmin <- 1
    pval_gene$xmax <- 2
    pval_gene <- pval_gene %>% add_significance(p.col = "padj")
  } 
  
  # Perform the t-test using the variable intgroup if no other groups are specified
  if ((is.null(select_group))){
    stat.test <- counts_data %>%
      t_test(as.formula(paste("count ~", intgroup)), ref.group = ref_group) %>%
      add_significance() %>%
      add_xy_position(x = intgroup)
  }
  
  # Reorder the factor levels for the group
  # counts_data[[intgroup]] <- factor(counts_data[[intgroup]], 
  #                                   levels = c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))
  # 
  # stat.test <- stat.test %>%
  #   mutate(
  #     xmin = as.numeric(factor(group1, levels = c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))),
  #     xmax = as.numeric(factor(group2, levels = c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))),
  #     y.position = c(max(counts_data) * 1.2, max(counts_data) * 1.4, max(counts_data) * 1.6)  # Adjust y positions dynamically
  #   )
  
  # Create a violin plot using ggplot2
  if ((is.null(select_group))){
    ggplot(counts_data, aes(x = .data[[intgroup]], y = count, fill = .data[[intgroup]])) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.6) +  # Adds individual points
      labs(title = paste(gene), y = "Normalized counts", x = "Group") +
      theme_classic() +
      scale_y_log10() +  # Optional: Log scale for better visualization +
      scale_fill_manual(values = c("blue", "red"))+
    theme(legend.position = "none") +
      stat_pvalue_manual(
        stat.test,
        label = "p.val.adj",
        tip.length = 0.01,label.size = 13,
        step.increase = 0.05,
        xmin = "xmin",
        xmax = "xmax",y.position = 1.1*log(max(counts_data$count),base = 10),
        inherit.aes = FALSE
      )+coord_cartesian(ylim =c(0.1*min(counts_data$count), 10*max(counts_data$count)))
  } else {
    ggplot(counts_data, aes(x = .data[[intgroup]], y = count, fill = .data[[intgroup]])) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.6) +  # Adds individual points
      labs(title = paste(gene), y = "Normalized counts", x = "Group") +
      theme_classic() +
      scale_y_log10() +  # Optional: Log scale for better visualization
      theme(legend.position = "none") +
      stat_pvalue_manual(
        pval_gene,label.size = 20,
        label = "padj.signif",
        tip.length = 0.01,
        step.increase = 0.05,
        xmin = "xmin",
        xmax = "xmax",y.position = 1.1*log(max(counts_data$count),base = 10),
        inherit.aes = FALSE
      )+coord_cartesian(ylim =c(0.1*min(counts_data$count), 10*max(counts_data$count)))
  }
}
