TCL_data_analysis
================
Mario Presti
First created on November 2024 Updated on 01 July 2026

- [Setup](#setup)
- [Dataset assembly](#dataset-assembly)
- [DE analysis and GSEA/Progeny](#de-analysis-and-gseaprogeny)
  - [DE analysis](#de-analysis)
  - [GSEA](#gsea)
  - [PROGENY and Dorothea](#progeny-and-dorothea)
- [Dedifferentiation and HLA-II
  positivity](#dedifferentiation-and-hla-ii-positivity)
  - [Tsoi states](#tsoi-states)
  - [Testing signature by Kim et al., JCI
    2021](#testing-signature-by-kim-et-al-jci-2021)
  - [Undifferentiated TCLs and HLA-II
    suppression](#undifferentiated-tcls-and-hla-ii-suppression)
  - [Inflammatory marker upregulation in HLA-II positive
    TCLs](#inflammatory-marker-upregulation-in-hla-ii-positive-tcls)
- [Differential HLA-II isotype
  expression](#differential-hla-ii-isotype-expression)
- [Save DDS](#save-dds)
- [Excluded analyses](#excluded-analyses)
  - [In-depth investigation of HLA-II phenotype in non-undifferentiated
    melanoma](#in-depth-investigation-of-hla-ii-phenotype-in-non-undifferentiated-melanoma)
  - [DE genes between HLA-II CONST+ and undifferentiated
    HLA-low](#de-genes-between-hla-ii-const-and-undifferentiated-hla-low)
  - [GSEA](#gsea-1)
  - [DE of specific inflammatory
    markers](#de-of-specific-inflammatory-markers)

# Setup

# Dataset assembly

This initial part of the script takes the transcriptomic data and
combines it into a single dataset

``` r
CCLE_counts <- read.csv(file = "../public/depmap_data/Depmap_melanoma.csv", header = T)
rownames(CCLE_counts) <- CCLE_counts$Features
rownames(CCLE_counts) <- sub("\\.","-",rownames(CCLE_counts))
CCLE_counts <- CCLE_counts[,-c(1,2)]
inHouse_counts <- read.table(file = "../controlled_data/TCL_data/salmon.merged.gene_counts.tsv", header = T)
#rename TMEM173 in STING1
rownames(inHouse_counts) <- inHouse_counts$gene_id
inHouse_counts["STING1",] <- inHouse_counts["TMEM173",]
inHouse_counts$gene_id <- NULL
inHouse_counts$gene_name <- NULL

Grasso_counts <- read.csv(file = "../public/grasso/d_raw_Grasso_RNAseq.csv")
rownames(Grasso_counts) <- Grasso_counts$X
Grasso_counts["STING1",] <- Grasso_counts["TMEM173",]
Grasso_counts$X <- NULL

CCLE_md <- read.csv(file = "../public/depmap_data/Depmap_md.csv", header = T)
rownames(CCLE_md) <- sub("-", "\\.",CCLE_md$ModelID)
CCLE_md <- CCLE_md[match(colnames(CCLE_counts), rownames(CCLE_md)), ]
CCLE_md$Dataset <- "CCLE"
CCLE_md$Sample <- rownames(CCLE_md)
inHouse_md <- read.xlsx(xlsxFile = "../controlled_data/TCL_data/Key_TCL.xlsx")
rownames(inHouse_md) <- inHouse_md$Sample
inHouse_md <- inHouse_md[match(colnames(inHouse_counts),rownames(inHouse_md)),]
inHouse_md$Dataset <- "inHouse"
inHouse_md$Sample <- rownames(inHouse_md)

#remove TCL15
inHouse_md <- subset(inHouse_md, Sample != "TCL15")

Grasso_md <- read.csv(file = "../public/grasso/md_raw_Grasso_RNAseq.csv")
rownames(Grasso_md) <- Grasso_md$Observation
Grasso_md <- Grasso_md[match(colnames(Grasso_counts), rownames(Grasso_md)),]
Grasso_md$Dataset <- "Grasso"
Grasso_md$Sample <- rownames(Grasso_md)

combined_md <- plyr::rbind.fill(CCLE_md, Grasso_md,inHouse_md)
rownames(combined_md) <- combined_md$Sample


CCLE_counts$gene <- rownames(CCLE_counts)
inHouse_counts$gene <- rownames(inHouse_counts)
Grasso_counts$gene <- rownames(Grasso_counts)
inHouse_counts <- inHouse_counts[,c(inHouse_md$Sample,"gene")]

combined_counts <- merge(CCLE_counts, Grasso_counts, by = "gene")
combined_counts <- merge(combined_counts, inHouse_counts, by="gene")

rm(list = c(grep("CCLE*",ls(.GlobalEnv), value = T),grep("Grasso*",ls(.GlobalEnv), value = T),grep("inHouse*",ls(.GlobalEnv), value = T)))
```

We perform unsupervised clustering of the samples based on their
expression of CIITA and HLA-DR genes. The first step is to visualize the
clustering output.

``` r
#heatmap pre-batch correction on normalized counts
rownames(combined_counts) <- combined_counts$gene
combined_counts$gene <- NULL

combined_dds <- DESeqDataSetFromMatrix(countData = round(combined_counts),
                                       colData = combined_md,
                                       design = ~ Dataset)
combined_vst <- vst(combined_dds, blind = TRUE)

# 2) Batch-correct the VST matrix for exploration/clustering
library(limma)
mat_vst_blind <- assay(combined_vst)
mat_vst_bc <- removeBatchEffect(mat_vst_blind, batch = combined_md$Dataset)

# 3) Wrap the corrected matrix back into a DESeqTransform-like object if you want
assay(combined_vst) <- mat_vst_bc

HLA_DR_genes <- grep("^HLA-DR", rownames(combined_dds), value = T)
HLA_DP_genes <- grep("^HLA-DP", rownames(combined_dds), value = T)
HLA_DQ_genes <- grep("^HLA-DQ", rownames(combined_dds), value = T)
HLAII_genes <- unique(c("CIITA", HLA_DR_genes, HLA_DP_genes, HLA_DQ_genes))
combined_clustering <- perform_cluster_analysis(combined_vst,HLAII_genes, annotation = "Dataset",perplexity = 30, clusters = 4)
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/clustering pre-batch correction-1.png" width="100%" />

``` r
combined_clustering$tSNE_plot
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/clustering pre-batch correction-2.png" width="100%" />

Dataset information+clustering was included in the DeSeq formula to
address the effects of the variables. In the end, we create a heatmap to
confirm that the clustering is isolating samples with high HLA
expression.The HLA-II genes have been removed from the heatmap.

``` r
#normalization, batch correction and clustering
#add the clustering info in the metadata
combined_md <- cbind(combined_md,combined_clustering$cluster_labels)
combined_md$KNN_cluster <- relevel(combined_md$KNN_cluster, ref = "HLA_low")
combined_dds <- DESeqDataSetFromMatrix(countData = round(combined_counts),
                                   colData = combined_md,
                                   design = ~ Dataset + KNN_cluster)
combined_dds <- DESeq(combined_dds)
#checking the differentially expressed genes
combined_results_filter <- as.data.frame(results(combined_dds, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(padj < 0.01 & baseMean >30 & abs(log2FoldChange) > 1)
combined_results_filter$gene <- rownames(combined_results_filter)
combined_genes <- rownames(combined_results_filter)

selected_ids <- c(combined_clustering$high_expression_ids, combined_clustering$low_expression_ids)

allgenes <- unique(rownames(combined_vst))
HLA_DR_genes <- grep("^HLA-DR", allgenes, value = T)
HLA_DP_genes <- grep("^HLA-DP", allgenes, value = T)
HLA_DQ_genes <- grep("^HLA-DQ", allgenes, value = T)
HLA_DO_genes <- grep("^HLA-DO", allgenes, value = T)
HLA_DM_genes <- grep("^HLA-DM", allgenes, value = T)
HLAII_genes <- unique(c("CIITA", "CD74", HLA_DR_genes,HLA_DP_genes,HLA_DQ_genes,HLA_DO_genes,HLA_DM_genes))
combined_genes <- allgenes[-which(allgenes %in% HLAII_genes)]
```

# DE analysis and GSEA/Progeny

## DE analysis

The differentially expressed genes can be visualized in a volcano plot,
where the HLA-II related genes have been removed, to focus on the other
transcriptomic differences between HLA-high and HLA-low samples

``` r
#preparing results for volcano plot without HLAII genes
combined_results_unfilter_no_HLA <- as.data.frame(results(combined_dds, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(abs(baseMean) >1000)
combined_results_unfilter_no_HLA$gene <- rownames(combined_results_unfilter_no_HLA)

combined_results_unfilter_no_HLA <- subset(combined_results_unfilter_no_HLA, !gene %in% HLAII_genes)
combined_results_unfilter_no_HLA <- subset(combined_results_unfilter_no_HLA,padj != "NA")
#volcano plot without HLA-II genes
ggplot(combined_results_unfilter_no_HLA, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, as.character(gene), ""),size = 3),
    size =5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  labs(title = "combined dataset - MHCII+ vs MHCII-",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/volcano plot-1.png" width="100%" />

``` r
#volcano plot with HLA-II genes
combined_results_unfilter <- as.data.frame(results(combined_dds, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(abs(baseMean) >300)
combined_results_unfilter$gene <- rownames(combined_results_unfilter)
ggplot(combined_results_unfilter, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, as.character(gene), ""),size = 3),
    size =5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  labs(title = "combined dataset - MHCII+ vs MHCII- including HLA-II genes",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/volcano plot-2.png" width="100%" />

``` r
diff_expr_genes <- combined_results_unfilter_no_HLA %>% filter(padj < 0.01 & abs(log2FoldChange) > 2 &  abs(baseMean) >300) %>%
  arrange(padj)
```

## GSEA

Then we perform GSEA to understand the phenotype of the HLA-II positive
cluster.

``` r
library(clusterProfiler)
library(org.Hs.eg.db)   
library(enrichplot)      

gene_list <- combined_results_unfilter_no_HLA$log2FoldChange
names(gene_list) <- combined_results_unfilter_no_HLA$gene
gene_list <- sort(gene_list, decreasing = TRUE)

msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
msigdb_hallmark_list <- msigdb_hallmark[, c("gs_name", "gene_symbol")]

# Run GSEA with the Hallmark gene sets
library(clusterProfiler)
gsea_results_plot <- GSEA(
  geneList = gene_list,           
  TERM2GENE = msigdb_hallmark_list,         
  pvalueCutoff = 0.5,            
  minGSSize = 10,                 
  maxGSSize = 500,                
  verbose = TRUE
)

# View the significant pathways with regulation status
library(ggrepel)
library(ggplot2)

# Convert NES to a bar plot friendly format by adding colors based on up- or down-regulation
gsea_results_plot@result$regulation <- ifelse(gsea_results_plot@result$NES > 0, "Upregulated", "Downregulated")

# Create the waterfall plot with bars
gsea_results_plot@result$ID <- sub("HALLMARK_", replacement = "", x = gsea_results_plot@result$ID)
ggplot(subset(gsea_results_plot@result,p.adjust<0.01), aes(x = reorder(ID, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = regulation), width = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  # Define the color for upregulated and downregulated pathways
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "GSEA Waterfall Plot - HLA+ vs HLA-",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)",
       fill = "Regulation") +
  theme(axis.text.y = element_text(size = 10, face = "bold"))  
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/GSEA-1.png" width="100%" />

``` r
ggplot(gsea_results_plot@result, aes(x = NES, y = -log10(p.adjust))) +
  geom_point(aes(color = ifelse(p.adjust < 0.01, "significant", "non-significant")), alpha = 0.6, size=3) +
  geom_text_repel(
    aes(label = ifelse(p.adjust < 0.05, as.character(ID), "")),
    size = 3,
    box.padding = unit(1.0, "lines"), 
    point.padding = unit(0, "lines"),
    segment.size = 0.5,  
    min.segment.length = 0.1, 
    nudge_y = 2,  
    nudge_x = 0.2  
  ) +
  theme_minimal() +
  labs(title = "GSEA combined Cell line dataset - HLA+ vs HLA-",
       x = "Normalized Enrichment Score",
       y = "-log10(p-value)",
       color = "Pathway Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/GSEA-2.png" width="100%" />

## PROGENY and Dorothea

``` r
library(progeny)
library(limma)
library(dplyr)

expr <- assay(combined_vst)
expr <- expr[!(rownames(expr) %in% HLAII_genes), , drop = FALSE]

path_act <- progeny(assay(combined_vst),
                    scale = TRUE,
                    organism = "Human",
                    top = 500,          
                    perm = 300,         
                    return_assay = TRUE)


# 2) Differential pathway activity: HLA_high vs HLA_low (adjust for dataset)
design <- model.matrix(~ Dataset + KNN_cluster, data = combined_md)
fit <- lmFit(t(path_act), design)            # pathways × samples -> transpose to samples × pathways
fit <- eBayes(fit)
progeny_res <- topTable(fit, coef = "KNN_clusterHLA_high", number = Inf) %>%
  tibble::rownames_to_column("Pathway")

sig <- progeny_res %>% dplyr::filter(adj.P.Val < 0.05)
# logFC > 0 => higher activity in HLA_high
ggplot(sig, aes(x = reorder(Pathway, logFC), y = logFC, fill = logFC > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "blue"), guide = "none") +
  labs(x = "Pathway", y = "Activity (logFC)", title = "PROGENy pathway activity: HLA+ vs HLA−")
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

The GSEA analysis shows that the HLA-II positive cell lines exhibit
signs of dedifferentiation, associated with upregulation of inflammatory
response markers

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/dotplot and GSEA plot-1.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/dotplot and GSEA plot-2.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/dotplot and GSEA plot-3.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/dotplot and GSEA plot-4.png" width="100%" />

# Dedifferentiation and HLA-II positivity

## Tsoi states

We also investigated if the Melanoma multi-stage differentiation
subtypes identified in Tsoi et al., Cancer Cell 2018 could be enriched
in HLA-II positive melanoma cell lines. Our analysis showed an
enrichment of the neural_crest-like phenotype, with a negative
enrichment for the Melanocytic phenotype, compatible with what has been
observed so far.

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/multi-stage differentiation subtypes-1.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/multi-stage differentiation subtypes-2.png" width="100%" />

## Testing signature by Kim et al., JCI 2021

As HLA-II follows dedifferentiation, how well is this recapitulated by
the signature in Kim et al, JCI 2021?
<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/Testing Kim et al signature-1.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/Testing Kim et al signature-2.png" width="100%" />

Finally, we show that the HLA-II CONST+ phenotype overlaps with the
Neural Crest-like phenotype

``` r
tsoi_map <- c(
  "Melanocytic"       = "Melan",
  "Transitory"        = "Transitory",
  "Neural_crest-like" = "Neur",
  "Undifferentiated"  = "Undiff"
)

tsoi_rows <- names(tsoi_map)

# 2) Subset Tsoi rows from ssGSEA matrix
tsoi_scores <- ssgsea_z[tsoi_rows, , drop = FALSE]

# 3) For each sample (column), pick the Tsoi state with the highest score
tsoi_call <- apply(tsoi_scores, 2, function(x) {
  tsoi_map[ tsoi_rows[ which.max(x) ] ]
})

# Quick sanity check: no NAs
table(is.na(tsoi_call))


# 4) Assign to metadata
combined_md <- combined_md %>% tibble::repair_names()
combined_md <- combined_md %>%
  dplyr::mutate(
    Tsoi = tsoi_call[Sample],
    Tsoi = factor(Tsoi, levels = c("Melan", "Transitory", "Neur", "Undiff"))
  )

table(combined_md$Tsoi, useNA = "ifany")


CIITA_counts <- assay(vst(combined_dds))["CIITA",]
melan_score <- ssgsea_z[ "Melanocytic",]
transit_score <- ssgsea_z[ "Transitory",]
neural_score <- ssgsea_z[ "Neural_crest-like",]
undiff_score <- ssgsea_z[ "Undifferentiated",]

melan_score <- melan_score[names(CIITA_counts)]
transit_score <- transit_score[names(CIITA_counts)]
neural_score <- neural_score[names(CIITA_counts)]
undiff_score <- undiff_score[names(CIITA_counts)]

scatterdata <- rbind(CIITA_counts, melan_score,transit_score,neural_score,undiff_score)
scatterdata <- as.data.frame(t(scatterdata))
scatterdata$Sample <- rownames(scatterdata)
scatterdata<- scatterdata %>% left_join(combined_md, by = "Sample")


mel<- ggscatter(x = "CIITA_counts",y =  "melan_score",
          add = "reg.line", conf.int = T,cor.coef = T,shape = "Dataset",
          scatterdata,title = "Melanocytic - CIITA expression" )
trans <- ggscatter(x = "CIITA_counts",y =  "transit_score",
          add = "reg.line", conf.int = T,cor.coef = T,shape = "Dataset",
          scatterdata, title = "Transitory - CIITA expression")
neur <- ggscatter(x = "CIITA_counts",y =  "neural_score",
          add = "reg.line", conf.int = T,cor.coef = T,shape = "Dataset",
          scatterdata, title = "Neural Crest-like - CIITA expression")
undiff <- ggscatter(x = "CIITA_counts",y =  "undiff_score",
          add = "reg.line", conf.int = T,cor.coef = T,shape = "Dataset",
          scatterdata, title = "Undifferentiated - CIITA expression")

ggarrange(mel,trans, neur,undiff)
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

## Undifferentiated TCLs and HLA-II suppression

In undifferentiated TCLs, IFNg is not able to upregulate HLA-II (Figure
S3)

``` r
cat(sub("\\.", "-",scatterdata[which(scatterdata$Tsoi == "Melan" & scatterdata$KNN_cluster == "HLA_low"),"Sample"]))
cat(sub("\\.", "-",scatterdata[which(scatterdata$Tsoi == "Undiff" & scatterdata$KNN_cluster == "HLA_low"),"Sample"]))
cat(sub("\\.", "-",scatterdata$Sample))

expr <- assay(combined_vst)
expr_mat <- as.matrix(expr)
mode(expr_mat) <- "numeric"

# (Optional) de-duplicate genes if needed (keep the most variable probe/gene)
if (any(duplicated(rownames(expr_mat)))) {
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
}

m_df <- msigdbr(species = "Homo sapiens", category = "H")

isg_df <- m_df %>%
  dplyr::filter(gs_name %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE"))

# Build list: 2 sets, each a vector of HGNC symbols
isg_list <- split(isg_df$gene_symbol, isg_df$gs_name)
isg_list$ISG_CORE <- c("JAK1", "JAK2", "STAT1", "IRF1", "IRF9")

### 2. Intersect with expressed genes (CRITICAL)
isg_list <- lapply(isg_list, function(g) base::intersect(g, rownames(expr_mat)))

sapply(isg_list, length)   # sanity check: both should be > 30–50 genes


## 3. Run ssGSEA to get per-sample ISG scores
# vst_mat should be a numeric matrix: rows = genes, cols = samples
gs <- gsvaParam(exprData = expr, 
          geneSets = isg_list)

isg_scores <- gsva(param = gs)
# isg_scores: matrix with rows = IFN gene sets, columns = samples


isg_df <- as.data.frame(t(isg_scores))
isg_df$Sample <- rownames(isg_df)

combined_md_new <- combined_md %>%
  left_join(isg_df, by = "Sample") %>%
  rename(
    IFNG_score = HALLMARK_INTERFERON_GAMMA_RESPONSE,
    IFNA_score = HALLMARK_INTERFERON_ALPHA_RESPONSE,
    INFLAMM_score = HALLMARK_INFLAMMATORY_RESPONSE
  )
## 5. Compare IFNγ scores across states

hla2_genes <- c("CIITA", "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74")
hla2_present <- base::intersect(hla2_genes, rownames(expr_mat))

hla2_mat   <- expr_mat[hla2_present, , drop = FALSE]
hla2_mat_z <- t(scale(t(hla2_mat)))

HLAII_score <- colMeans(hla2_mat_z, na.rm = TRUE)

combined_md_new$HLAII_score <- HLAII_score

combined_md_new %>%
  dplyr::group_by(Tsoi) %>%
  dplyr::summarise(
    n          = n(),
    cor_IFN_HL2 = cor(INFLAMM_score, HLAII_score, use = "complete.obs")
  )

combined_md_new$Tsoi <- factor(combined_md$Tsoi,
                           levels = c("Melan","Transitory","Neur","Undiff"))

fit <- lm(HLAII_score ~ INFLAMM_score * Tsoi, data = combined_md_new)
summary(fit)

combined_md_new$Tsoi <- relevel(combined_md_new$Tsoi, "Neur")
fit <- lm(HLAII_score ~ INFLAMM_score * Tsoi, data = combined_md_new)
summary(fit)



ggplot(combined_md_new, aes(x = IFNG_score, y = HLAII_score, colour = Tsoi)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(
    method      = "pearson",
    size        = 5, label.y = 1.5
  ) +
  facet_wrap(~Tsoi) +
  labs(
    x = "IFN/ISG score",
    y = "HLA-II module score",
    title = "State-dependent coupling between IFN activity and HLA-II"
  ) +
  theme_bw() +
  scale_color_manual(
    values = c(
      Melan      = "#0073C2FF",
      Transitory = "#EFC000FF",
      Neur       = "#CD534CFF",
      Undiff     = "#4A6990FF"
    )
  )
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/State-dependent coupling between IFN activity and HLA-II-1.png" width="100%" />

## Inflammatory marker upregulation in HLA-II positive TCLs

And we show potential mechanistic drivers relying in increased
inflammatory potential (exploratory)

``` r
progeny_df <- as.data.frame(path_act)
progeny_df$Sample <- rownames(progeny_df)

# join to your existing scatterdata
scatterdata2 <- scatterdata %>%
  left_join(progeny_df, by = "Sample")

path_cols <- c("NFkB", "Androgen", "EGFR", "Estrogen", "JAK.STAT",
               "MAPK", "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF",
               "WNT")  # adapt to your names


cors <- sapply(path_cols, function(p) {
  cor(scatterdata2[[p]], scatterdata2$CIITA_counts, use = "complete.obs", method = "pearson")
})

cor_df <- data.frame(
  pathway = path_cols,
  r       = cors
)

scatter_nc <- scatterdata2 %>% dplyr::filter(neural_score > median(neural_score))  # or > median

cors_nc <- sapply(path_cols, function(p) {
  cor(scatter_nc[[p]], scatter_nc$CIITA_counts, use = "complete.obs")
})
cor_nc_df <- data.frame(pathway = path_cols, r = cors_nc)


library(ggplot2)

ggplot(cor_nc_df, aes(x = reorder(pathway, r), y = r)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "",
    y = "Pearson R with CIITA",
    title = "Pathway activity correlated with CIITA in NC-like melanomas"
  ) +
  theme_bw()
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/potential inflammatory drivers p-1.png" width="100%" />

``` r
fit <- lm(
  CIITA_counts ~ neural_score + NFkB + TNFa + JAK.STAT,
  data = scatterdata2
)
summary(fit)

scatter_scaled <- scatterdata2 %>%
  mutate(
    CIITA_z      = scale(CIITA_counts)[,1],
    neural_z     = scale(neural_score)[,1],
    NFkB_z       = scale(NFkB)[,1],
    TNFa_z       = scale(TNFa)[,1],
    JAKSTAT_z    = scale(JAK.STAT)[,1]
  )

fit_z <- lm(CIITA_z ~ neural_z + NFkB_z + TNFa_z + JAKSTAT_z, data = scatter_scaled)
summary(fit_z)$coefficients

coef_df <- broom::tidy(fit_z, conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(
    term = factor(
      term,
      levels = c("neural_z", "JAKSTAT_z", "NFkB_z", "TNFa_z"),
      labels = c("Neural crest score", "JAK–STAT activity",
                 "NFκB activity", "TNFα activity")
    )
  )

p_coef <- ggplot(coef_df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  labs(
    x = "Standardized effect on CIITA (β ± 95% CI)",
    y = "",
    title = "Relative contribution of melanoma state and inflammatory pathways to CIITA"
  ) +
  theme_bw()

p_coef
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/potential inflammatory drivers p-2.png" width="100%" />

# Differential HLA-II isotype expression

We also investigated the differential expression of HLA-II isotypes in
our cohort and it showed a specific a higher expression of HLA-DR
compared to the other isotypes

``` r
genes_to_plot <- c("HLA-DRA", "HLA-DPA1", "HLA-DQA1")
combined_dds$KNN_cluster <- factor(combined_dds$KNN_cluster, levels = c("HLA_low", "HLA_mid_low", "HLA_mid_high", "HLA_high"))
create_violin_plot_compare_genes_vst(
  dds = combined_dds,
  genes = genes_to_plot,
  intgroup = "KNN_cluster",
)+scale_fill_manual(values = c("#0072B2", "#E69F00","#009E40"))+theme(axis.text=element_text(size=12),
                                                                      axis.title=element_text(size=14,face="bold"),
                                                                      axis.text.x = element_text(size = 15, face = "bold"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/HLA isotype comparison-1.png" width="100%" />

(Not included) Code used to generate the genes used for string analysis

``` r
cat(subset(combined_results_filter, log2FoldChange >= 2 & padj < 0.01 & baseMean >100)[,"gene"])
```

# Save DDS

``` r
save(combined_dds, file="../controlled_data/combined_dds.RData")
save(combined_clustering, file="../controlled_data/combined_clustering.RData")
```

# Excluded analyses

## In-depth investigation of HLA-II phenotype in non-undifferentiated melanoma

Differences between HLA-II CONST+ and melanocytic HLA-low

``` r
# undiff_tcls <- names(which(ssgsea_scores["Undifferentiated",] > ssgsea_group_z["Undifferentiated",]))
undiff_hla_low <- which(ids_undiff %in% combined_clustering$low_expression_ids)
undiff_hla_high <- which(ids_undiff %in% combined_clustering$high_expression_ids)

rate_low <- length(undiff_hla_low)/length(combined_clustering$low_expression_ids)
rate_high <- length(undiff_hla_high)/length(combined_clustering$high_expression_ids)

combined_md$undifferentiated <- ifelse(combined_md$Sample %in% ids_undiff, "Yes", "No")

low_md <- subset(combined_md, KNN_cluster == "HLA_low" & undifferentiated == "No")
low_md <- rbind(low_md,subset(combined_md, KNN_cluster != "HLA_low"))

low_counts <- combined_counts[,low_md$Sample]

undiff_dds <- DESeqDataSetFromMatrix(countData = round(combined_counts),
                                       colData = combined_md,
                                       design = ~ Dataset + KNN_cluster)
undiff_dds <- DESeq(undiff_dds)
#checking "all" the differentially expressed genes
undiff_results_filter <- as.data.frame(results(undiff_dds, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(padj < 0.01 & baseMean >30 & abs(log2FoldChange) > 1)
undiff_results_filter$gene <- rownames(undiff_results_filter)

ggplot(undiff_results_filter, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, as.character(gene), ""),size = 3),
    size =5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  labs(title = "HLA low - Undifferentiated vs Melanocytic",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/analyses without undiff melanoma to exclude confounding due to diff status-1.png" width="100%" />

``` r
undiff_results_filter_no_HLA <- as.data.frame(results(undiff_dds, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(abs(baseMean) >300)
undiff_results_filter_no_HLA$gene <- rownames(undiff_results_filter_no_HLA)

undiff_results_filter_no_HLA <- subset(undiff_results_filter_no_HLA, !gene %in% HLAII_genes)
undiff_results_filter_no_HLA <- subset(undiff_results_filter_no_HLA,padj != "NA")
#volcano plot without HLA-II genes
ggplot(undiff_results_filter_no_HLA, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, as.character(gene), ""),size = 3),
    size =5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  labs(title = "combined dataset - MHCII+ vs MHCII- melanocytic subset",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/analyses without undiff melanoma to exclude confounding due to diff status-2.png" width="100%" />

## DE genes between HLA-II CONST+ and undifferentiated HLA-low

``` r
# undiff_tcls <- names(which(ssgsea_scores["Undifferentiated",] > ssgsea_group_z["Undifferentiated",]))
low_md_und <- subset(combined_md, KNN_cluster == "HLA_low" & undifferentiated == "Yes")
low_md_und <- rbind(low_md_und,subset(combined_md, KNN_cluster != "HLA_low"))

low_counts_und <- combined_counts[,low_md_und$Sample]

undiff_dds_undiff <- DESeqDataSetFromMatrix(countData = round(low_counts_und),
                                       colData = low_md_und,
                                       design = ~ Dataset + KNN_cluster)
undiff_dds_undiff <- DESeq(undiff_dds_undiff)
#checking "all" the differentially expressed genes
undiff_results_und_filter <- as.data.frame(results(undiff_dds_undiff, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(padj < 0.01 & baseMean >30 & abs(log2FoldChange) > 1)
undiff_results_und_filter$gene <- rownames(undiff_results_und_filter)

ggplot(undiff_results_und_filter, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, as.character(gene), ""),size = 3),
    size =5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  labs(title = "MHC-IICONST+ vs MHCII-CONST- Undifferentiated subset",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/DE HLA-II const+ and undiff HLA-low-1.png" width="100%" />

``` r
undiff_results_und_filter_no_HLA <- as.data.frame(results(undiff_dds_undiff, name = "KNN_cluster_HLA_high_vs_HLA_low",cooksCutoff = T)) %>% 
  dplyr::filter(abs(baseMean) >300)
undiff_results_und_filter_no_HLA$gene <- rownames(undiff_results_und_filter_no_HLA)

undiff_results_und_filter_no_HLA <- subset(undiff_results_und_filter_no_HLA, !gene %in% HLAII_genes)
undiff_results_und_filter_no_HLA <- subset(undiff_results_und_filter_no_HLA,padj != "NA")
#volcano plot without HLA-II genes
ggplot(undiff_results_und_filter_no_HLA, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(padj < 0.01 & abs(log2FoldChange) > 2, as.character(gene), ""),size = 3),
    size =5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  labs(title = "combined dataset - MHC-IICONST+ vs MHCII-CONST- Undifferentiated subset",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/DE HLA-II const+ and undiff HLA-low-2.png" width="100%" />

## GSEA

``` r
gene_list_melan <- undiff_results_filter_no_HLA$log2FoldChange
names(gene_list_melan) <- undiff_results_filter_no_HLA$gene
gene_list_melan <- sort(gene_list_melan, decreasing = TRUE)

# Run GSEA with the Hallmark gene sets
gsea_results_melan <- GSEA(
  geneList = gene_list_melan,           # Ranked gene list
  TERM2GENE = msigdb_hallmark_list,          # Hallmark gene sets in the correct format
  pvalueCutoff = 0.5,            # P-value cutoff for significance
  minGSSize = 10,                 # Minimum gene set size
  maxGSSize = 500,                # Maximum gene set size
  verbose = TRUE
)

# Assuming gsea_results contains the GSEA results with columns NES, p.adjust, and ID

# Convert NES to a bar plot friendly format by adding colors based on up- or down-regulation
gsea_results_melan@result$regulation <- ifelse(gsea_results_melan@result$NES > 0, "Upregulated", "Downregulated")

# Create the waterfall plot with bars
gsea_results_melan@result$ID <- sub("HALLMARK_", replacement = "", x = gsea_results_melan@result$ID)
ggplot(subset(gsea_results_melan@result,p.adjust<0.05), aes(x = reorder(ID, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = regulation), width = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  coord_flip() + 
  theme_minimal() +
  labs(title = "GSEA Waterfall Plot - HLA-IICONST+ vs Melanocytic HLA-IICONST-",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)",
       fill = "Regulation") +
  theme(axis.text.y = element_text(size = 10, face = "bold"))  
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/GSEA on undiff dds-1.png" width="100%" />

``` r
ggplot(gsea_results_melan@result, aes(x = NES, y = -log10(p.adjust))) +
  geom_point(aes(color = ifelse(p.adjust < 0.05, "significant", "non-significant")), alpha = 0.6, size=3) +
  geom_text_repel(
    aes(label = ifelse(p.adjust < 0.05, as.character(ID), "")),
    size = 3,
    box.padding = unit(1.0, "lines"), 
    point.padding = unit(0, "lines"),
    segment.size = 0.5,  
    min.segment.length = 0.1,force = 3,
    nudge_y = 2,  
    nudge_x = 0.2  
  ) +
  theme_minimal() +
  labs(title = "GSEA combined Cell line dataset - HLA-IICONST+ vs Melanocytic HLA-IICONST-",
       x = "Normalized Enrichment Score",
       y = "-log10(p-value)",
       color = "Pathway Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/GSEA on undiff dds-2.png" width="100%" />

``` r
gene_list_undiff <- undiff_results_und_filter_no_HLA$log2FoldChange
names(gene_list_undiff) <- undiff_results_und_filter_no_HLA$gene
gene_list_undiff <- sort(gene_list_undiff, decreasing = TRUE)

# Run GSEA with the Hallmark gene sets
gsea_results_undiff <- GSEA(
  geneList = gene_list_undiff,           # Ranked gene list
  TERM2GENE = msigdb_hallmark_list,          # Hallmark gene sets in the correct format
  pvalueCutoff = 0.5,            # P-value cutoff for significance
  minGSSize = 10,                 # Minimum gene set size
  maxGSSize = 500,                # Maximum gene set size
  verbose = TRUE
)

# Convert NES to a bar plot friendly format by adding colors based on up- or down-regulation
gsea_results_undiff@result$regulation <- ifelse(gsea_results_undiff@result$NES > 0, "Upregulated", "Downregulated")

# Create the waterfall plot with bars
gsea_results_undiff@result$ID <- sub("HALLMARK_", replacement = "", x = gsea_results_undiff@result$ID)
ggplot(subset(gsea_results_undiff@result,p.adjust<0.05), aes(x = reorder(ID, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = regulation), width = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  coord_flip() +  
  theme_minimal() +
  labs(title = "GSEA Waterfall Plot - HLA-IICONST+ vs Undifferentiated HLA-IICONST-",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)",
       fill = "Regulation") +
  theme(axis.text.y = element_text(size = 10, face = "bold"))  
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/GSEA on undiff dds-3.png" width="100%" />

``` r
ggplot(gsea_results_undiff@result, aes(x = NES, y = -log10(p.adjust))) +
  geom_point(aes(color = ifelse(p.adjust < 0.05, "significant", "non-significant")), alpha = 0.6, size=3) +
  geom_text_repel(
    aes(label = ifelse(p.adjust < 0.05, as.character(ID), "")),
    size = 3,
    box.padding = unit(1.0, "lines"), 
    point.padding = unit(0, "lines"),
    segment.size = 0.5,  
    min.segment.length = 0.1, force = 3,
    nudge_y = 2,  
    nudge_x = 0.2  
  ) +
  theme_minimal() +
  labs(title = "GSEA combined Cell line dataset - HLA-IICONST+ vs Undifferentiated HLA-IICONST-",
       x = "Normalized Enrichment Score",
       y = "-log10(p-value)",
       color = "Pathway Significance") +
  scale_color_manual(values = c("grey", "red"))
```

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/GSEA on undiff dds-4.png" width="100%" />

## DE of specific inflammatory markers

Finally, to confirm these results we studied the differences in
expression of specific genes, linked with both dedifferentiation and
inflammatory responses. Our results show that HLA-II positive melanomas
have a constitutive activation of NFKB1, probably driving both
dedifferentiation and CIITA activation

<img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-1.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-2.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-3.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-4.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-5.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-6.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-7.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-8.png" width="100%" /><img src="E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/TCL_analysis_final_files/figure-gfm/violin plots on undiff - combined dds-9.png" width="100%" />
