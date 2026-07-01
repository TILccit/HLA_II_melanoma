---
editor_options: 
  markdown: 
    wrap: 72
---

# Genomic_drivers_HLA_II

Mario Presti First created on November 2024 Updated on 01 July 2026

-   [Genomic differences between HLA+ and HLA-
    TCLs](#genomic-differences-between-hla-and-hla--tcls)
-   [TMB differences in HLA+ and HLA-
    TCLs](#tmb-differences-in-hla-and-hla--tcls)
-   [Checking the distribution of NRAS and BRAF
    muts](#checking-the-distribution-of-nras-and-braf-muts)
-   [Genomic signatures related to HLA-II
    positivity](#genomic-signatures-related-to-hla-ii-positivity)
-   [Duplication of specific genes and HLA
    positivity](#duplication-of-specific-genes-and-hla-positivity)

# Genomic differences between HLA+ and HLA- TCLs {#genomic-differences-between-hla-and-hla--tcls}

``` r
#trying to identify a mutational driver of CIITA activation

library(maftools)
coding_classes <- c(
  "Missense_Mutation",  "Nonsense_Mutation",
  "Frame_Shift_Del",     "Frame_Shift_Ins",
  "In_Frame_Del",       "In_Frame_Ins"
)
mutations <- read.maf(maf = "public/depmap_data/OmicsSomaticMutationsMAFProfile.maf",vc_nonSyn = coding_classes)
```

```         
## -Reading
## -Validating
## --Removed 4 duplicated variants
## -Silent variants: 38896 
## -Summarizing
## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
##   AHNAK2
##   OBSCN
## -Processing clinical data
## --Missing clinical data
## -Finished in 9.840s elapsed (8.890s cpu)
```

``` r
mutations_HLAsamples_high <- subsetMaf(mutations, tsb = c(profile_id_wgs_high$ProfileID))
```

```         
## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
## -Processing clinical data
```

``` r
mutations_HLAsamples_low <- subsetMaf(mutations, tsb = c(profile_id_wgs_low$ProfileID))
```

```         
## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
## -Processing clinical data
```

``` r
compar <- mafCompare(m1 = mutations_HLAsamples_high, m2 = mutations_HLAsamples_low, m1Name = 'HLA+', m2Name = 'HLA-', minMut = 5)
```

# TMB differences in HLA+ and HLA- TCLs {#tmb-differences-in-hla-and-hla--tcls}

``` r
capture_mb <- 30
tmb_high <- tmb(mutations_HLAsamples_high, captureSize = 30)
```

```         
## Filtering CNV events (if any..)

## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
## -Processing clinical data
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
tmb_low  <- tmb(mutations_HLAsamples_low, captureSize = 30)
```

```         
## Filtering CNV events (if any..)

## --Possible FLAGS among top ten genes:
##   TTN
##   MUC16
## -Processing clinical data
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# 8. compare distributions
library(ggplot2)
df <- rbind(
  data.table(group = "HLA+", tmb_high),
  data.table(group = "HLA-",  tmb_low)
)

ggplot(df, aes(x = group, y = total_perMB_log)) +
  geom_boxplot() + geom_point()+
  ylab("log₂(TMB per Mb)") +
  theme_minimal()
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

# Checking the distribution of NRAS and BRAF muts {#checking-the-distribution-of-nras-and-braf-muts}

##NRAS

``` r
gene_of_interest <- "NRAS"  

# 2) get a flat data.table of all samples in the MAF
all_samples <- data.table(
  Tumor_Sample_Barcode = getSampleSummary(mutations)$Tumor_Sample_Barcode
)

maf_dt <- data.table::as.data.table(mutations@data)
mut_flag <- maf_dt[
  Hugo_Symbol == gene_of_interest, 
  .(HasMut = TRUE), 
  by = Tumor_Sample_Barcode
]

status_dt <- merge(
  all_samples, 
  mut_flag, 
  by = "Tumor_Sample_Barcode", 
  all.x = TRUE
)
status_dt[, HasMut := ifelse(is.na(HasMut), FALSE, TRUE)]

hla_dt <- data.table(
  Tumor_Sample_Barcode = c(profile_id_wgs_high$ProfileID, profile_id_wgs_low$ProfileID),
  HLA_status = c(
    rep("High", length(profile_id_wgs_high$ProfileID)),
    rep("Low",  length(profile_id_wgs_low$ProfileID))
  )
)
plot_dt_NRAS <- merge(status_dt, hla_dt, by = "Tumor_Sample_Barcode")

# 5) prepare a factor for plotting
plot_dt_NRAS[, MutStatus := factor(
  ifelse(HasMut, "Mutated", "Wild‑type"),
  levels = c("Wild‑type","Mutated")
)]

rate_dt <- plot_dt_NRAS[
  , .(MutRate = mean(HasMut)),    # fraction mutated
  by = HLA_status
]

diff_rate_NRAS <- rate_dt[HLA_status=="High", MutRate] -
  rate_dt[HLA_status=="Low",  MutRate]

rate_dt_NRAS <- rate_dt

rate_dt_NRAS$HLA_status <- as.factor(rate_dt_NRAS$HLA_status)
rate_dt_NRAS$HLA_status <- relevel(rate_dt_NRAS$HLA_status, ref="Low")


ggplot(rate_dt_NRAS, aes(x = HLA_status, y = MutRate, fill = HLA_status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = percent(MutRate)),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(), limits = c(0,1)) +
  labs(
    x = "HLA-II status",
    y = "Mutation rate",
    title = paste0(gene_of_interest, " mutation rate by HLA status"),
    subtitle = paste0("Δ rate (High – Low) = ", percent(diff_rate_NRAS)),
    caption = paste0("N high=", length(which(rate_dt_NRAS$HLA_status=="High")), "; N low=", length(which(rate_dt_NRAS$HLA_status=="Low")))
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

##BRAF

``` r
gene_of_interest <- "BRAF"  

# 2) get a flat data.table of all samples in the MAF
all_samples <- data.table(
  Tumor_Sample_Barcode = getSampleSummary(mutations)$Tumor_Sample_Barcode
)

maf_dt <- data.table::as.data.table(mutations@data)
mut_flag <- maf_dt[
  Hugo_Symbol == gene_of_interest, 
  .(HasMut = TRUE), 
  by = Tumor_Sample_Barcode
]

status_dt <- merge(
  all_samples, 
  mut_flag, 
  by = "Tumor_Sample_Barcode", 
  all.x = TRUE
)
status_dt[, HasMut := ifelse(is.na(HasMut), FALSE, TRUE)]

hla_dt <- data.table(
  Tumor_Sample_Barcode = c(profile_id_wgs_high$ProfileID, profile_id_wgs_low$ProfileID),
  HLA_status = c(
    rep("High", length(profile_id_wgs_high$ProfileID)),
    rep("Low",  length(profile_id_wgs_low$ProfileID))
  )
)
plot_dt <- merge(status_dt, hla_dt, by = "Tumor_Sample_Barcode")

# 5) prepare a factor for plotting
plot_dt[, MutStatus := factor(
  ifelse(HasMut, "Mutated", "Wild‑type"),
  levels = c("Wild‑type","Mutated")
)]

rate_dt <- plot_dt[
  , .(MutRate = mean(HasMut)),    # fraction mutated
  by = HLA_status
]

diff_rate <- rate_dt[HLA_status=="High", MutRate] -
  rate_dt[HLA_status=="Low",  MutRate]

rate_dt$HLA_status <- as.factor(rate_dt$HLA_status)
rate_dt$HLA_status <- relevel(rate_dt$HLA_status, ref="Low")

ggplot(rate_dt, aes(x = HLA_status, y = MutRate, fill = HLA_status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = percent(MutRate)),
            vjust = -0.5, size = 5) +
  scale_y_continuous(labels = percent_format(), limits = c(0,1)) +
  labs(
    x = "HLA-II status",
    y = "Mutation rate",
    title = paste0(gene_of_interest, " mutation rate by HLA status"),
    subtitle = paste0("Δ rate (High – Low) = ", percent(diff_rate)),
    caption = paste0("N high=", length(which(plot_dt$HLA_status=="High")), "; N low=", length(which(plot_dt$HLA_status=="Low")))
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Genomic signatures related to HLA-II positivity {#genomic-signatures-related-to-hla-ii-positivity}

``` r
signatures <- fread(paste0("public/depmap_data/OmicsSignatures.csv"))
colnames(signatures)[1] <- "ModelID"
signatures$ModelID <- sub("-", ".", signatures$ModelID)
signatures_HLA <- subset(signatures, ModelID %in% c(model_used_high, model_used_low))
signatures_HLA$HLA_status <- ifelse(signatures_HLA$ModelID %in% model_used_high,"HLA_high", "HLA_low")
signatures_HLA$HLA_status <- as.factor(signatures_HLA$HLA_status) %>% relevel(ref = "HLA_low")
ggplot(signatures_HLA, aes(x = HLA_status, y = CIN, fill = HLA_status)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(
    title = "Violin plot of CIN by HLA status",
    y = "CIN",
    x = "HLA-II status"
  ) +
  theme_minimal() +
  scale_y_log10() +
  scale_fill_manual(values = c("blue", "red")) +
  theme(legend.position = "none")+ 
   stat_compare_means(
  method      = "t.test",
  label       = "p.format",
  comparisons = list(c("HLA_low", "HLA_high")),
  label.y.npc = "top",
  size        = 5
)
```

```         
## Warning: Removed 11 rows containing non-finite outside the scale range
## (`stat_ydensity()`).

## Warning: Removed 11 rows containing non-finite outside the scale range
## (`stat_signif()`).

## Warning: Removed 11 rows containing missing values or values outside the scale range
## (`geom_point()`).
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggplot(signatures_HLA, aes(x = HLA_status, y = MSIScore, fill = HLA_status)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(
    title = "Violin plot of MSI score by HLA status",
    y = "MSIScore",
    x = "HLA-II status"
  ) +
  theme_minimal() +
  scale_y_log10() +
  scale_fill_manual(values = c("blue", "red")) +
  theme(legend.position = "none")+
  stat_compare_means(
  method      = "t.test",
  label       = "p.format",
  comparisons = list(c("HLA_low", "HLA_high")),
  label.y.npc = "top",
  size        = 5
)
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggplot(signatures_HLA, aes(x = HLA_status, y = LoHFraction, fill = HLA_status)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(
    title = "Violin plot of LoHFraction score by HLA status",
    y = "LoHFraction",
    x = "HLA-II status"
  ) +
  theme_minimal() +
  scale_y_log10() +
  scale_fill_manual(values = c("blue", "red")) +
  theme(legend.position = "none")+
     stat_compare_means(
  method      = "t.test",
  label       = "p.format",
  comparisons = list(c("HLA_low", "HLA_high")),
  label.y.npc = "top",
  size        = 5
)
```

```         
## Warning: Removed 11 rows containing non-finite outside the scale range
## (`stat_ydensity()`).

## Warning: Removed 11 rows containing non-finite outside the scale range
## (`stat_signif()`).

## Warning: Removed 11 rows containing missing values or values outside the scale range
## (`geom_point()`).
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

# Duplication of specific genes and HLA positivity {#duplication-of-specific-genes-and-hla-positivity}

``` r
#could it be linked to duplication of specific genes?
CNV_data <- fread("public/depmap_data/OmicsCNGene.csv")
colnames(CNV_data)[1] <- "ModelID"
CNV_data$ModelID <- sub("-", ".", CNV_data$ModelID)
my_models <- c(model_used_high,model_used_low)
CNV_data_sub <- CNV_data[ModelID %chin% my_models] 
# %chin% is a faster version of %in% for character vectors

gene_of_interest <- "CIITA"
i <- grep(gene_of_interest,colnames(CNV_data_sub))[1]
gene_CNV <- CNV_data_sub[,c(1,..i)]

# Add factor for HLA-high/low
gene_CNV[, HLA_class := ifelse(ModelID %in% model_used_high, "high", "low")]
col_name <- sym(colnames(gene_CNV)[grep(gene_of_interest, colnames(gene_CNV))])
# Make the boxplot
ggplot(gene_CNV, aes(x = HLA_class, y= !!col_name, fill = HLA_class)) +
  geom_boxplot()
```

![](E:/PhD_projects/MHCII_project/Paper/Clean_scripts/GitHub/reports/Genomic_drivers_HLA_II_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
