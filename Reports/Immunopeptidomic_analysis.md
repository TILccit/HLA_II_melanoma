Immunopeptidomics_analysis
================
Mario Presti
First created on November 2024 Updated on 01 July 2026

- [Setup](#setup)
- [Data parsing and QC](#data-parsing-and-qc)
- [IFNg induced peptides](#ifng-induced-peptides)
- [Exclusivity scores](#exclusivity-scores)
  - [Calculation and plotting for untreated
    conditions](#calculation-and-plotting-for-untreated-conditions)
  - [Calculation and plotting for IFNg treated
    conditions](#calculation-and-plotting-for-ifng-treated-conditions)
- [NetMHCIIPan analysis](#netmhciipan-analysis)
  - [HLA-II isotype redundancy
    analysis](#hla-ii-isotype-redundancy-analysis)
  - [UMAP projection of NetMHCIIpan](#umap-projection-of-netmhciipan)
  - [Score distribution in each
    cluster](#score-distribution-in-each-cluster)
- [CAPtan analysis](#captan-analysis)
- [Captan results](#captan-results)
- [Core motif analysis (Exploratory, not
  included)](#core-motif-analysis-exploratory-not-included)

# Setup

# Data parsing and QC

This initial part of the script takes the transcriptomic data and
combines it into a single dataset

``` r
tbl <- as.data.frame(with(res, table(Master.Gene.name)))

ggplot(as.data.frame(tbl), aes(x=Master.Gene.name, y=Freq))+geom_col(position = 'dodge')+
  geom_label_repel(data=subset(tbl, Freq > 100), aes(label = Master.Gene.name) )
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/Data%20parsing%20and%20cleaning-1.png)<!-- -->

``` r
ggplot(tbl, aes(x = Master.Gene.name, y = Freq)) +
  geom_col(fill = "steelblue", width = 5) +
  geom_text(
    aes(label = Master.Gene.name),
    vjust = -0.5,                      # push text just above the bar
    size = 3,                          # adjust text size if needed
    angle = 45,                        # rotate if labels collide
    hjust = 0,
    data = subset(tbl, Freq > 100)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),    # hide crowded x-axis labels
    axis.ticks.x = element_blank()
  )
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/Data%20parsing%20and%20cleaning-2.png)<!-- -->

``` r
wt_col <- grep("WT_",colnames(res))
wt_ifn_col <- grep("WT..IFN",colnames(res))

dr_col <- grep("DR.KO_R",colnames(res))
dr_ifn_col <- grep("DR.KO..IFN",colnames(res))

dp_col <- grep("DP.KO_R",colnames(res))
dp_ifn_col <- grep("DP.KO..IFN",colnames(res))

dq_col <- grep("DQ.KO_R",colnames(res))
dq_ifn_col <- grep("DQ.KO..IFN",colnames(res))

ciita_col <- grep("CIITA",colnames(res))

res$AVG_untreated_WT <- apply(res[,wt_col], MARGIN = 1, mean)
res$AVG_IFN_WT <- apply(res[,wt_ifn_col], MARGIN = 1, mean)

#filter out all the peptides expressed  more in the CIITA condition
res$AVG_CIITA <- apply(res[,ciita_col], MARGIN = 1, mean)

library(forcats)
library(tidyr)
#Compare WT to controls

#pull out ordering unique sequence for peptides
res$order <- order(res$AVG_untreated_WT)


res_long <- res %>%
  dplyr::select(order, AVG_untreated_WT, AVG_IFN_WT, AVG_CIITA) %>%
  pivot_longer(
    cols = c(AVG_untreated_WT, AVG_IFN_WT, AVG_CIITA),
    names_to = "Condition",
    values_to = "Expression"
  ) %>%
  mutate(
    Condition = recode(Condition,
                       AVG_untreated_WT = "WT",
                       AVG_IFN_WT = "WT_IFN",
                       AVG_CIITA        = "CIITA"
    )
  )


WT_vs_ctrls <- ggplot(res_long, 
                      aes(
                        x=order,
                        y = Expression,
                        fill = Condition)) +
  geom_col(width = 5, alpha=0.5, position= "identity") +
  theme_minimal() +
  scale_fill_manual(
    name   = "Condition",
    values = c("WT"    = "lightblue",
               "WT_IFN" = "tomato",
               "CIITA" = "black")
  ) +
  theme_bw(base_size = 14) +                        # white background
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"                         # move legend if you like
  ) +
  labs(
    x = NULL,
    y = "Average fold enrichment"
  )

#correlation plots
library(smplot2)
```

    ## Updated tutorial for smplot2: smin95.github.io/dataviz/

``` r
WT_vs_IFN_corr <- ggplot(res, (aes(x = AVG_untreated_WT, y = AVG_IFN_WT))) +
  geom_point() +
  sm_statCorr()+
  theme_minimal(base_size = 14) +
  labs(
    x = "Average WT",
    y = "Average WT IFN"
  )
```

    ## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## ℹ The deprecated feature was likely used in the smplot2 package.
    ##   Please report the issue at <https://github.com/smin95/smplot2/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
WT_vs_CIITA_corr <- ggplot(res, (aes(x = AVG_untreated_WT, y = AVG_CIITA))) +
  geom_point() +
  sm_statCorr()+
  theme_minimal(base_size = 14) +
  labs(
    x = "Average WT",
    y = "Average CIITA"
  )

#label the peptides where CIITA > WT_IFN
res$exclude <- ifelse(
  (res$AVG_CIITA - res$AVG_untreated_WT) > 0,
  "QC pass","QC not pass"
)

WT_vs_IFN_corr_split <- ggplot(res, (aes(x = AVG_untreated_WT, y = AVG_IFN_WT))) +
  geom_point() +
  facet_grid(cols = vars(exclude))+
  sm_statCorr()+
  theme_minimal(base_size = 14) +
  labs(
    x = "Average WT",
    y = "Average WT IFN"
  )

#### subsetting only QC passed ####
res <- subset(res, exclude == "QC pass")
```

# IFNg induced peptides

What are the differences between IFNg and noIFNg conditions in terms of
isolated peptides?

``` r
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(limma)
})

# --- 0) Select WT and IFN columns (WT_… vs WT…IFN per your regex) ---
wt_cols  <- grep("^WT", names(res), value = TRUE)
ifn_cols <- grep("^WT.*IFN",  names(res), value = TRUE)
wt_cols <- setdiff(wt_cols, ifn_cols)
stopifnot(length(wt_cols) > 0, length(ifn_cols) > 0)

E <- as.matrix(res[, c(wt_cols, ifn_cols)])
mode(E) <- "numeric"
rownames(E) <- res$Sequence
gene_of <- res$Master.Gene.name

## 2) Summarize to gene x sample by robust median over peptides for that gene
samples <- colnames(E)
G <- sapply(samples, function(s)
  tapply(E[, s], INDEX = gene_of, FUN = function(v) median(v, na.rm = TRUE))
)
G <- as.matrix(G); mode(G) <- "numeric"

## 3) Filter genes with too few finite values in either group
wt_idx  <- match(wt_cols,  colnames(G))
ifn_idx <- match(ifn_cols, colnames(G))
keep <- rowSums(is.finite(G[, wt_idx,  drop=FALSE])) >= 2 &
  rowSums(is.finite(G[, ifn_idx, drop=FALSE])) >= 2
G <- G[keep, , drop=FALSE]

## 4) limma design and fit (IFN − UT)
group  <- factor(c(rep("UT", length(wt_cols)), rep("IFN", length(ifn_cols))),
                 levels = c("UT","IFN"))
group <- relevel(group,ref = "IFN")
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit   <- lmFit(G, design)                     # continuous data; NA allowed
fitc  <- contrasts.fit(fit, makeContrasts(IFN - UT, levels = design))
fitc  <- eBayes(fitc, robust = TRUE, trend = TRUE)

gene_tt <- topTable(fitc, number = Inf, sort.by = "P")
# gene_tt$logFC > 0 → higher in IFN; < 0 → higher in untreated
gene_tt$gene <- rownames(gene_tt)
# #plot as volcano
ggplot(gene_tt, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1.5, "significant", "non-significant")), alpha = 0.6) +
  geom_text_repel(
    aes(label = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1.5, as.character(gene), ""),size = 3),
    size =5,max.overlaps = 25,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() + ylim(0,3)+
  labs(title = "combined dataset - MHCII+ vs MHCII-",
       x = "Log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Significance") +
  scale_color_manual(values = c("grey", "red"))
```

    ## Warning: ggrepel: 92 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/ifn%20-%20no%20ifn%20differences-1.png)<!-- -->

# Exclusivity scores

## Calculation and plotting for untreated conditions

``` r
res$DR <- apply(res[,dr_col], MARGIN = 1, mean)
res$DR_IFN <- apply(res[,dr_ifn_col], MARGIN = 1, mean)

#calculate DP score
res$DP <- apply(res[,dp_col], MARGIN = 1, mean)
res$DP_IFN <- apply(res[,dp_ifn_col], MARGIN = 1, mean)

#calculate DQ score
res$DQ <- apply(res[,dq_col], MARGIN = 1, mean)
res$DQ_IFN <- apply(res[,dq_ifn_col], MARGIN = 1, mean)

res <- res %>% 
  dplyr::mutate(
    # DR metrics
    DR_depletion   = (abs(DR   - AVG_CIITA)),
    DR_exclusivity = ((DP  - DR) +
                        (DQ   - DR))   / 2,
    DR_score       = DR_exclusivity*4 - DR_depletion,
    
    # DP metrics
    DP_depletion   = (abs(DP   - AVG_CIITA)),
    DP_exclusivity = ((DR   - DP) +
                        (DQ   - DP))   / 2,
    DP_score       = DP_exclusivity*4 - DP_depletion,
    
    # DQ metrics
    DQ_depletion   = (abs(DQ   - AVG_CIITA)),
    DQ_exclusivity = ((DR   - DQ) +
                        (DP   - DQ))   / 2,
    DQ_score       = DQ_exclusivity*4 - DQ_depletion
  )

pep_order <- res %>% 
  group_by(Annotated.Sequence) %>% 
  dplyr::summarize(med_DR = median(DR_score, na.rm=TRUE)) %>% 
  arrange(desc(med_DR)) %>% 
  pull(Annotated.Sequence)

library(tibble)
library(pheatmap)
# 2) pivot your scores to long form - choose the row with the best scores
score_mat <- res |>
  dplyr::select(Sequence, DR_score, DP_score, DQ_score) |>
  dplyr::group_by(Sequence) |>
  dplyr::slice_max(order_by = pmax(abs(DR_score), abs(DP_score), abs(DQ_score)),
                   n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  tibble::column_to_rownames("Sequence") |>
  as.matrix()

ggplot(as.data.frame(score_mat),aes(x=DR_score,y=DP_score))+
  geom_point()+geom_smooth(method="glm")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20untreated-1.png)<!-- -->

``` r
ggplot(as.data.frame(score_mat),aes(x=DR_score,y=DQ_score))+
  geom_point()+geom_smooth(method="glm")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20untreated-2.png)<!-- -->

``` r
ggplot(as.data.frame(score_mat),aes(x=DP_score,y=DQ_score))+
  geom_point()+geom_smooth(method="glm")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20untreated-3.png)<!-- -->

``` r
# score_mat <- ramify::clip(score_mat,.min = -2)
# 2) Draw the heatmap with clustering
pheatmap(
  score_mat,
  scale         = "row",     # standardize each row (peptide) to mean=0, sd=1
  clustering_distance_rows  = "euclidean",
  clustering_method         = "ward.D2",
  cluster_cols  = FALSE,     # if you only want to cluster peptides, not alleles
  show_rownames = FALSE,     # hide peptide names if there are hundreds
  color         = colorRampPalette(c("steelblue", "white", "tomato"))(100),
  main          = "Private-peptide scores (clustered)"
)
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20untreated-4.png)<!-- -->

``` r
#use k-means clustering to divide the peptides into classes
set.seed(1234)         # for reproducibility
k <- 3               # number of clusters
km <- kmeans(score_mat, centers = k, nstart = 25)

# km$cluster gives a vector of length = nrow(score_mat)
res <- res %>%
  mutate(
    KM_Cluster = km$cluster[Sequence]
  )

# Quick table
table(res$KM_Cluster)
```

    ## 
    ##    1    2    3 
    ## 3567 6777  511

``` r
#now plot them in a umap
library(uwot)
```

    ## Loading required package: Matrix
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
# 2) Run UMAP on the 3‐dimensional score space
set.seed(1234)
umap_emb <- umap(
  score_mat,
  n_neighbors = 15,
  min_dist     = 0.1,
  metric       = "euclidean"
)

# 3) Build a data.frame of UMAP coords + cluster labels
df_umap <- as.data.frame(umap_emb)
colnames(df_umap) <- c("UMAP1", "UMAP2")
df_umap$Peptide <- rownames(score_mat)

# pull in your k-means cluster per peptide
cluster_map <- res %>%
  dplyr::distinct(Sequence, KM_Cluster) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(cluster_map, by = "Peptide")

#add the HLA-DR score
DR_score_map <- res %>%
  distinct(Sequence, DR_score) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(DR_score_map, by = "Peptide")

#add the HLA-DP score
DP_score_map <- res %>%
  distinct(Sequence, DP_score) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(DP_score_map, by = "Peptide")
```

    ## Warning in left_join(., DP_score_map, by = "Peptide"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 28 of `x` matches multiple rows in `y`.
    ## ℹ Row 2132 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
#add the HLA-DQ score
DQ_score_map <- res %>%
  distinct(Sequence, DQ_score) %>%
  dplyr::rename(Peptide = Sequence)

df_umap_no_ifn <- df_umap %>%
  left_join(DQ_score_map, by = "Peptide")
```

    ## Warning in left_join(., DQ_score_map, by = "Peptide"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 28 of `x` matches multiple rows in `y`.
    ## ℹ Row 2132 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
# 4) Plot, coloring by KM_Cluster
umap_no_ifn <- ggplot(df_umap_no_ifn, aes(x = UMAP1, y = UMAP2, color = factor(KM_Cluster))) +
  geom_point(alpha = 0.8, size = 1.5) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "UMAP of peptide privateness scores"
  ) +
  scale_color_viridis_d(option = "C")

DR_umap_noifn <- ggplot(df_umap_no_ifn, aes(x = UMAP1, y = UMAP2, color = DR_score)) +
  geom_point(alpha = 0.8, size = 1.5,show.legend = F) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",,
    title  = "HLA-DR score"
  ) +
  scale_color_viridis_c(option = "inferno")


DP_umap_noifn <- ggplot(df_umap_no_ifn, aes(x = UMAP1, y = UMAP2, color = DP_score)) +
  geom_point(alpha = 0.8, size = 1.5,show.legend = F) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "HLA-DP score"
  ) +
  scale_color_viridis_c(option = "inferno")

DQ_umap_noifn <- ggplot(df_umap_no_ifn, aes(x = UMAP1, y = UMAP2, color = DQ_score)) +
  geom_point(alpha = 0.8, size = 1.5, show.legend = F) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "HLA-DQ score"
  ) +
  scale_color_viridis_c(option = "inferno")

library(ggpubr)

plotlist <- list(DR_umap_noifn, DP_umap_noifn, DQ_umap_noifn)
ggarrange(plotlist = plotlist, common.legend = T, legend = "bottom")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20untreated-5.png)<!-- -->

## Calculation and plotting for IFNg treated conditions

``` r
res <- res %>% 
  mutate(
    # DR metrics
    DR_depletion   = (abs(DR_IFN   - AVG_CIITA)),
    DR_exclusivity = ((DP_IFN   - DR_IFN) +
                        (DQ_IFN   - DR_IFN))   / 2,
    DR_score       = DR_exclusivity*4 - DR_depletion,
    
    # DP metrics
    DP_depletion   = (abs(DP_IFN   - AVG_CIITA)),
    DP_exclusivity = ((DR_IFN   - DP_IFN) +
                        (DQ_IFN   - DP_IFN))   / 2,
    DP_score       = DP_exclusivity*4 - DP_depletion,
    
    # DQ metrics
    DQ_depletion   = (abs(DQ_IFN   - AVG_CIITA)),
    DQ_exclusivity = ((DR_IFN   - DQ_IFN) +
                        (DP_IFN   - DQ_IFN))   / 2,
    DQ_score       = DQ_exclusivity*4 - DQ_depletion
  )

pep_order_ifn <- res %>% 
  group_by(Annotated.Sequence) %>% 
  dplyr::summarize(med_DR = median(DR_score, na.rm=TRUE)) %>% 
  arrange(desc(med_DR)) %>% 
  pull(Annotated.Sequence)

library(tibble)
library(pheatmap)
# 2) pivot your scores to long form - choose the row with the best scores
score_mat_ifn <- res |>
  dplyr::select(Sequence, DR_score, DP_score, DQ_score) |>
  dplyr::group_by(Sequence) |>
  dplyr::slice_max(order_by = pmax(abs(DR_score), abs(DP_score), abs(DQ_score)),
                   n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  tibble::column_to_rownames("Sequence") |>
  as.matrix()

ggplot(as.data.frame(score_mat_ifn),aes(x=DR_score,y=DP_score))+
  geom_point()+geom_smooth(method="glm")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20treated-1.png)<!-- -->

``` r
ggplot(as.data.frame(score_mat_ifn),aes(x=DR_score,y=DQ_score))+
  geom_point()+geom_smooth(method="glm")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20treated-2.png)<!-- -->

``` r
ggplot(as.data.frame(score_mat_ifn),aes(x=DP_score,y=DQ_score))+
  geom_point()+geom_smooth(method="glm")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20treated-3.png)<!-- -->

``` r
#use k-means clustering to divide the peptides into classes
set.seed(1234)         # for reproducibility
k <- 3               # number of clusters
km <- kmeans(score_mat_ifn, centers = k, nstart = 25)

# km$cluster gives a vector of length = nrow(score_mat)
res <- res %>%
  mutate(
    KM_Cluster = km$cluster[Sequence]
  )

# Quick table
table(res$KM_Cluster)
```

    ## 
    ##    1    2    3 
    ## 2001 6460 2394

``` r
#now plot them in a umap
library(uwot)


# 2) Run UMAP on the 3‐dimensional score space
set.seed(1234)
umap_emb <- umap(
  score_mat_ifn,
  n_neighbors = 15,
  min_dist     = 0.1,
  metric       = "euclidean"
)

# 3) Build a data.frame of UMAP coords + cluster labels
df_umap <- as.data.frame(umap_emb)
colnames(df_umap) <- c("UMAP1", "UMAP2")
df_umap$Peptide <- rownames(score_mat_ifn)

# pull in your k-means cluster per peptide
cluster_map <- res %>%
  dplyr::distinct(Sequence, KM_Cluster) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(cluster_map, by = "Peptide")

#add the HLA-DR score
DR_score_map <- res %>%
  distinct(Sequence, DR_score) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(DR_score_map, by = "Peptide")

#add the HLA-DP score
DP_score_map <- res %>%
  distinct(Sequence, DP_score) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(DP_score_map, by = "Peptide")
```

    ## Warning in left_join(., DP_score_map, by = "Peptide"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 28 of `x` matches multiple rows in `y`.
    ## ℹ Row 2132 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
#add the HLA-DQ score
DQ_score_map <- res %>%
  distinct(Sequence, DQ_score) %>%
  dplyr::rename(Peptide = Sequence)

df_umap <- df_umap %>%
  left_join(DQ_score_map, by = "Peptide")
```

    ## Warning in left_join(., DQ_score_map, by = "Peptide"): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 28 of `x` matches multiple rows in `y`.
    ## ℹ Row 2132 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
# 4) Plot, coloring by KM_Cluster
ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = factor(KM_Cluster))) +
  geom_point(alpha = 0.8, size = 1.5) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "UMAP of peptide privateness scores"
  ) +
  scale_color_viridis_d(option = "C")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20treated-4.png)<!-- -->

``` r
DR_umap <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = DR_score)) +
  geom_point(alpha = 0.8, size = 1.5,show.legend = F) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",,
    title  = "HLA-DR score"
  ) +
  scale_color_viridis_c(option = "inferno")


DP_umap <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = DP_score)) +
  geom_point(alpha = 0.8, size = 1.5,show.legend = F) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "HLA-DP score"
  ) +
  scale_color_viridis_c(option = "inferno")

DQ_umap <- ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = DQ_score)) +
  geom_point(alpha = 0.8, size = 1.5, show.legend = F) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "HLA-DQ score"
  ) +
  scale_color_viridis_c(option = "inferno")

plotlist_ifn <- list(DR_umap, DP_umap, DQ_umap)
ggarrange(plotlist = plotlist_ifn, common.legend = T, legend = "bottom")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/exclusivity%20treated-5.png)<!-- -->

``` r
#write results
res_umap <- merge(res,df_umap,by.x = "Sequence", by.y = "Peptide")

#write results to put them in netMHCIIpan
write.table(res_umap$Sequence, file = "../controlled_data/Immunopeptidomics/all.peptides.tsv", append = F,quote = F)
```

# NetMHCIIPan analysis

First, we parse the NetMHCIIPan outputs

``` r
netMHCIIres <- parse_netmhciipan_xls("../controlled_data/Immunopeptidomics/netMHCII_out/all_rank2.xls", strong_thr = 1, weak_thr = 5)
netMHCIIres<- subset(netMHCIIres, Binder != "NB")

netMHCIIres$Gene <- ifelse(grepl("DR", netMHCIIres$Allele), "HLA_DR", 
                           ifelse(grepl("DP", netMHCIIres$Allele), "HLA_DP",
                                  ifelse(grepl("DQ", netMHCIIres$Allele), "HLA_DQ",NA)))
netMHCIIres_tbl <- as.data.frame(prop.table(table(netMHCIIres$Gene, netMHCIIres$Binder)))
colnames(netMHCIIres_tbl) <- c("Gene", "Binder", "Proportion")
netMHCIIres_tbl$Gene <- factor(netMHCIIres_tbl$Gene, levels = c("HLA_DR", "HLA_DP", "HLA_DQ"))

ggplot(netMHCIIres_tbl, aes(x=Gene))+
  geom_col(aes(y=Proportion, fill=Binder), position = "dodge2")+theme_classic()
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/Parse%20NetMHCIIPan%20output-1.png)<!-- -->

## HLA-II isotype redundancy analysis

Second, we test the redundancy of each peptide in the different alleles

``` r
weak_thr <- 5  # set 1 for SB-only analysis

redundancy <- netMHCIIres |>
  dplyr::filter(Rank <= weak_thr) |>
  dplyr::distinct(Peptide, Allele) |>
  dplyr::count(Peptide, name = "n_alleles_binding") |>
  dplyr::arrange(desc(n_alleles_binding))

# quick view of redundancy distribution
redundancy |> dplyr::count(n_alleles_binding)
```

    ## # A tibble: 6 × 2
    ##   n_alleles_binding     n
    ##               <int> <int>
    ## 1                 1  3809
    ## 2                 2  1220
    ## 3                 3   445
    ## 4                 4   144
    ## 5                 5   103
    ## 6                 6    47

``` r
# how is this redundancy distributed?
library(tidyr)
x_gene <- netMHCIIres |>
  mutate(Gene = dplyr::case_when(
    grepl("^DR", Allele) ~ "HLA_DR",
    grepl("^HLA-DP", Allele) ~ "HLA_DP",
    grepl("^HLA-DQ", Allele) ~ "HLA_DQ",
    TRUE ~ NA_character_
  ))

# binary flags per peptide × gene family (binder if EL rank <= weak_thr)
bind_flags <- x_gene |>
  dplyr::filter(Rank <= weak_thr) |>
  distinct(Peptide, Gene) |>
  mutate(bind = 1L) |>
  tidyr::pivot_wider(names_from = Gene, values_from = bind, values_fill = 0)

# counts of binding patterns across families
pattern_counts <- bind_flags |>
  dplyr::mutate(pattern = paste0(HLA_DR, HLA_DP, HLA_DQ)) |>
  dplyr::count(pattern, name = "n_peptides") |>
  dplyr::arrange(desc(n_peptides))
pattern_counts
```

    ## # A tibble: 7 × 2
    ##   pattern n_peptides
    ##   <chr>        <int>
    ## 1 100           2367
    ## 2 010           1899
    ## 3 101            594
    ## 4 110            318
    ## 5 001            268
    ## 6 111            214
    ## 7 011            108

``` r
#Among peptides whose best allele is DP or DQ, how many also bind DR (≤ weak_thr)?
  
best_call <- x_gene |>
  group_by(Peptide) |>
  slice_min(Rank, with_ties = FALSE) |>    
  ungroup() |>
  transmute(Peptide, best_Allele = Allele, best_Gene = Gene, best_Rank = Rank)
  
dpdq_best <- best_call |> dplyr::filter(best_Gene %in% c("HLA_DP", "HLA_DQ"))
  
dpdq_with_DR <- x_gene |>
   dplyr::filter(Rank <= weak_thr, Gene == "HLA_DR") |>
   distinct(Peptide) |>
   inner_join(dpdq_best, by = "Peptide")
  
frac_dpdq_also_DR <- nrow(dpdq_with_DR) / nrow(dpdq_best)
frac_dpdq_also_DR
```

    ## [1] 0.1469816

``` r
# For all DR binders, what fraction also bind any DP-DQ?
DR_best <- best_call |> dplyr::filter(best_Gene %in% c("HLA_DR"))

DR_binders <- x_gene |>
  dplyr::filter(Rank <= weak_thr, Gene %in% c("HLA_DR")) |>
   distinct(Peptide) |>
  inner_join(DR_best, by = "Peptide")
  
DR_binders_with_dpdq <- DR_binders |>
   inner_join(x_gene |> dplyr::filter(Rank <= weak_thr, Gene %in% c("HLA_DP", "HLA_DQ")) |> distinct(Peptide),
              by = "Peptide")
  
frac_any_DR_also_dpdq <- nrow(DR_binders_with_dpdq) / nrow(DR_binders)
frac_any_DR_also_dpdq
```

    ## [1] 0.2366978

## UMAP projection of NetMHCIIpan

``` r
per_gene_minrank <- x_gene |>
  group_by(Peptide, Gene) |>
  dplyr::summarize(minRank = min(Rank, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Gene, values_from = minRank,
                     names_prefix = "minRank_", values_fill = NA)

binder_flags <- x_gene |>
  dplyr::filter(Rank <= 5) |>
  distinct(Peptide, Gene) |>
  mutate(flag = 1L) |>
  tidyr::pivot_wider(names_from = Gene, values_from = flag,
                     names_prefix = "bind_", values_fill = 0)

netmhc_summary <- per_gene_minrank |>
  full_join(binder_flags, by = "Peptide")

# harmonize the key name if needed and full outer join
merged <- res_umap |>
  dplyr::rename(Peptide = Sequence) |>
  dplyr::full_join(netmhc_summary, by = "Peptide")

# if you prefer to keep the original column name 'peptide'
merged <- merged |>
  dplyr::rename(peptide = Peptide)
merged_sub <- merged[which(!is.na(merged$bind_HLA_DP)),]

ggplot(merged,aes(x = DR_score.y, y = minRank_HLA_DR))+geom_point()+geom_smooth(method = "glm", se = FALSE)
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 11102 rows containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning: Removed 11102 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/NetMHCIIpan%20UMAPs-1.png)<!-- -->

``` r
merged$bind_HLA_DP <- as.factor(merged$bind_HLA_DP)
merged$bind_HLA_DQ <- as.factor(merged$bind_HLA_DQ)
merged$bind_HLA_DR <- as.factor(merged$bind_HLA_DR)

ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = factor(KM_Cluster))) +
  geom_point(alpha = 0.8, size = 1.5) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "UMAP of peptide privateness scores"
  ) +
  scale_color_viridis_d(option = "C")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/NetMHCIIpan%20UMAPs-2.png)<!-- -->

``` r
ggplot(merged, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.8, size = 1.5,aes(color=bind_HLA_DR) ) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "UMAP of peptide privateness scores"
  )
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/NetMHCIIpan%20UMAPs-3.png)<!-- -->

``` r
ggplot(merged, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.8, size = 1.5,aes(color=bind_HLA_DP) ) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "UMAP of peptide privateness scores"
  )
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/NetMHCIIpan%20UMAPs-4.png)<!-- -->

``` r
ggplot(merged, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.8, size = 1.5,aes(color=bind_HLA_DQ) ) +
  theme_bw(base_size = 14) +
  labs(
    x      = "UMAP 1",
    y      = "UMAP 2",
    color  = "k-means Cluster",
    title  = "UMAP of peptide privateness scores"
  )
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/NetMHCIIpan%20UMAPs-5.png)<!-- -->

## Score distribution in each cluster

``` r
cluster_dr <- subset(res, KM_Cluster == 3 & DQ_score < quantile(res$DQ_score)[4])
cluster_dp <- subset(res, KM_Cluster == 1 & DQ_score < quantile(res$DQ_score)[4])
cluster_dq <- subset(res, KM_Cluster %in% c(2) & DQ_score > quantile(res$DQ_score)[3] & DR_score < quantile(res$DR_score,0.5) & DP_score < quantile(res$DP_score,0.5))

res <- res %>% mutate(final_clusters =
                        case_when(
                          res$KM_Cluster == 3 & DQ_score < quantile(res$DQ_score)[4] ~ "DR",
                          res$KM_Cluster == 1 & DQ_score < quantile(res$DQ_score)[4]  ~ "DP",
                          res$KM_Cluster %in% c(2) & DQ_score > quantile(res$DQ_score)[3] & DR_score < quantile(res$DR_score,0.5) & DP_score < quantile(res$DP_score,0.5) ~ "DQ",
                          TRUE ~ NA_character_
                          ))
#plot distribution of scores in clusters
res_long <- res %>%
  dplyr::select(Sequence, final_clusters, DR_score, DP_score, DQ_score) %>%
  pivot_longer(cols = c(DR_score, DP_score, DQ_score),
               names_to = "Score_type",
               values_to = "Value")

# Drop NAs if needed
res_long <- res_long %>% dplyr::filter(!is.na(final_clusters))
res_long$final_clusters <- factor(res_long$final_clusters, levels=c("DR","DP","DQ")) 
res_long$Score_type <- factor(res_long$Score_type, levels=c("DR_score","DP_score","DQ_score")) 

res_merged <- res %>%
  left_join(merged %>% dplyr::select(peptide, bind_HLA_DR, bind_HLA_DP, bind_HLA_DQ),
            by = c("Sequence" = "peptide"))
```

    ## Warning in left_join(., merged %>% dplyr::select(peptide, bind_HLA_DR, bind_HLA_DP, : Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 18 of `x` matches multiple rows in `y`.
    ## ℹ Row 1305 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
# define binder assignment
bind_labels <- merged %>%
  mutate(
    DR_bind = ifelse(!is.na(bind_HLA_DR) & bind_HLA_DR == "1", 1, 0),
    DP_bind = ifelse(!is.na(bind_HLA_DP) & bind_HLA_DP == "1", 1, 0),
    DQ_bind = ifelse(!is.na(bind_HLA_DQ) & bind_HLA_DQ == "1", 1, 0)
  ) %>%
  rowwise() %>%
  mutate(Binder_label = {
    total <- DR_bind + DP_bind + DQ_bind
    if (total == 1) {
      if (DR_bind == 1) "DR"
      else if (DP_bind == 1) "DP"
      else "DQ"
    } else if (total == 0) {
      "NA"   # nonbinder
    } else {
      NA_character_  # multiple conflicting binders
    }
  }) %>%
  ungroup() %>%
  dplyr::select(peptide, Binder_label)

res_long <- res_long %>%
  left_join(bind_labels, by = c("Sequence" = "peptide"))
```

    ## Warning in left_join(., bind_labels, by = c(Sequence = "peptide")): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 25 of `x` matches multiple rows in `y`.
    ## ℹ Row 5756 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

``` r
table(res_long$Binder_label, res_long$final_clusters)
```

    ##     
    ##         DR    DP    DQ
    ##   DP    30  9303   114
    ##   DQ   348     9  1425
    ##   DR 10038   513   264
    ##   NA  3276  1023  1947

``` r
#rename to the original cluster name
res_long$final_clusters <- ifelse(res_long$final_clusters == "DR", 3,
                                  ifelse(res_long$final_clusters == "DP", 1,2))

# Plot distributions (violin + boxplot)
ggplot(res_long, aes(x = Score_type, y = Value, fill = Score_type)) +
  geom_point(aes(color=Binder_label), alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  facet_wrap(~ final_clusters, scales = "free_y") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  xlab("Cluster (DR, DP, DQ)") +
  ylab("Score value")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/Score-cluster%20histograms-1.png)<!-- -->

``` r
binder_summary <- res_long %>%
  dplyr::filter(!is.na(final_clusters)) %>%
  # dplyr::filter(Binder_label != "NA") %>%
  group_by(final_clusters, Binder_label) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  group_by(final_clusters) %>%
  mutate(prop = n / sum(n))

# Stacked proportion barplot
ggplot(binder_summary, aes(x = final_clusters, y = prop, fill = Binder_label)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = "bottom"
  ) +
  ylab("Proportion of predicted binders") +
  xlab("Cluster (DR, DP, DQ)") 
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/Score-cluster%20histograms-2.png)<!-- -->

# CAPtan analysis

Export results for captan

``` r
write_fasta_with_genes(
  sequences = cluster_dp$Sequence,
  genes = cluster_dp$Master.Gene.name,
  file = "../controlled_data/Immunopeptidomics/captan/DP_peptides_with_genes.fasta"
)

write_fasta_with_genes(
  sequences = cluster_dr$Sequence,
  genes = cluster_dr$Master.Gene.name,
  file = "../controlled_data/Immunopeptidomics/captan/DR_peptides_with_genes.fasta"
)

write_fasta_with_genes(
  sequences = cluster_dq$Sequence,
  genes = cluster_dq$Master.Gene.name,
  file = "../controlled_data/Immunopeptidomics/captan/DQ_peptides_with_genes.fasta"
)
```

# Captan results

``` r
dr_cap <- read.csv("../controlled_data/Immunopeptidomics/captan/out_DR/summary_CAPTAn_ctx.csv", header = T)
dp_cap <- read.csv("../controlled_data/Immunopeptidomics/captan/out_DP/summary_CAPTAn_ctx.csv", header = T)
dq_cap <- read.csv("../controlled_data/Immunopeptidomics/captan/out_DQ/summary_CAPTAn_ctx.csv", header = T)
dr_cap$gene <- "HLA-DR"
dp_cap$gene <- "HLA-DP"
dq_cap$gene <- "HLA-DQ"
cap_all <- rbind(dr_cap, dp_cap, dq_cap)
cap_all$gene <- factor(cap_all$gene, levels=c("HLA-DR", "HLA-DP", "HLA-DQ")) 

thresh<-quantile(cap_all$Score)[3]
ggplot(cap_all, aes(x=gene, y=Score, fill=gene)) +
  geom_violin(trim=FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, outlier.shape=NA, color="black") +
  geom_jitter(width=0.1, alpha=0.05, size=0.5) +
  geom_hline(yintercept=thresh, linetype="dashed") +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("HLA-DR", "HLA-DP"), c("HLA-DR", "HLA-DQ")),
                     label = "p.signif",size=14)+
  theme(legend.position="none", axis.text.x = element_text(size = 14)) +
  ylab("CAPTAn antigenicity score") + xlab("")+scale_y_continuous(expand = 0.1)
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
dr_top <- subset(dr_cap, Score > quantile(dr_cap$Score)[4])
dp_top <- subset(dp_cap, Score > quantile(dp_cap$Score)[4])
dq_top <- subset(dq_cap, Score > quantile(dq_cap$Score)[4])
top_all <- rbind(dr_top, dp_top,dq_top)
```

# Core motif analysis (Exploratory, not included)

``` r
core_dr <- read.csv(file = "../controlled_data/Immunopeptidomics/captan/out_DR/summary_CAPTAn_core.csv")
core_dr_filt_first <- subset(core_dr[which(core_dr$Protein %in% dr_top$Protein),], Primary_core != "" & Score > quantile(core_dr$Score,0.9))
ggseqlogo(core_dr_filt_first$Primary_core, seq_type="aa") +
  theme_bw(base_size=14) +
  labs(title="Logo of central 9-mers from top DR peptides")
```

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the ggseqlogo package.
    ##   Please report the issue at <https://github.com/omarwagih/ggseqlogo/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/motif%20analysis-1.png)<!-- -->

``` r
core_dr_filt_second <- subset(core_dr[which(core_dr$Protein %in% dr_top$Protein),], Secondary_core != "" & Score > quantile(core_dr$Score,0.9))
ggseqlogo(core_dr_filt_second$Secondary_core, seq_type="aa") +
  theme_bw(base_size=14) +
  labs(title="Logo of central 9-mers from top DR peptides (second core)")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/motif%20analysis-2.png)<!-- -->

``` r
core_dp <- read.csv(file = "../controlled_data/Immunopeptidomics/captan/out_DP/summary_CAPTAn_core.csv")
core_dp_filt_first <- subset(core_dp[which(core_dp$Protein %in% dp_top$Protein),], Primary_core != "" & Score > quantile(core_dp$Score,0.9))
ggseqlogo(core_dp_filt_first$Primary_core, seq_type="aa") +
  theme_bw(base_size=14) +
  labs(title="Logo of central 9-mers from top DP peptides")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/motif%20analysis-3.png)<!-- -->

``` r
core_dp_filt_second <- subset(core_dp[which(core_dp$Protein %in% dp_top$Protein),], Secondary_core != "" & Score > quantile(core_dp$Score,0.9))
ggseqlogo(core_dp_filt_second$Secondary_core, seq_type="aa") +
  theme_bw(base_size=14) +
  labs(title="Logo of central 9-mers from top DP peptides (second core)")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/motif%20analysis-4.png)<!-- -->

``` r
core_dq <- read.csv(file = "../controlled_data/Immunopeptidomics/captan/out_DQ/summary_CAPTAn_core.csv")
core_dq_filt_first <- subset(core_dq, Primary_core != "" & Score > quantile(core_dq$Score,0.9))
ggseqlogo(core_dq_filt_first$Primary_core, seq_type="aa") +
  theme_bw(base_size=14) +
  labs(title="Logo of central 9-mers from top DQ peptides")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/motif%20analysis-5.png)<!-- -->

``` r
core_dq_filt_second <- subset(core_dq, Secondary_core != "" & Score > quantile(core_dq$Score,0.9))
ggseqlogo(core_dq_filt_second$Secondary_core, seq_type="aa") +
  theme_bw(base_size=14) +
  labs(title="Logo of central 9-mers from top DQ peptides (second core)")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/Reports/Immunopeptidomic_analysis_files/figure-gfm/motif%20analysis-6.png)<!-- -->
