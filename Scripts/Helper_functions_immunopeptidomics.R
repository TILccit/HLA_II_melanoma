parse_netmhciipan_xls <- function(path, strong_thr = 1, weak_thr = 5) {
  # Read everything as character (safer for messy "xls")
  raw <- read_tsv(path, col_names = FALSE,
                  col_types = cols(.default = col_character()),
                  progress = FALSE)
  
  # Headers
  h1 <- as.character(raw[1, ]); h2 <- as.character(raw[2, ])
  
  # Remove columns that are completely empty (header+body)
  keep <- vapply(seq_len(ncol(raw)), function(j) {
    a <- h1[j]; b <- h2[j]
    if (!is.na(a) && nzchar(a)) return(TRUE)
    if (!is.na(b) && nzchar(b)) return(TRUE)
    vals <- raw[-c(1,2), j, drop = TRUE]
    any(!(is.na(vals) | vals == ""))
  }, logical(1))
  raw <- raw[, keep, drop = FALSE]
  h1 <- as.character(raw[1, ]); h2 <- as.character(raw[2, ])
  
  # Forward-fill allele names across blocks
  grp <- h1; grp[grp==""] <- NA
  last <- NA_character_
  for (j in seq_along(grp)) {
    if (!is.na(grp[j])) last <- grp[j]
    grp[j] <- last
  }
  
  # Impute missing metric names within each allele block
  metrics_order <- c("Core","Inverted","Score","Rank","Score_BA","nM","Rank_BA")
  idx_within <- integer(length(grp))
  cur <- NA_character_; k <- 0L
  for (j in seq_along(grp)) {
    g <- grp[j]
    if (is.na(g)) { idx_within[j] <- NA_integer_; next }
    if (!identical(g, cur)) { cur <- g; k <- 0L }
    k <- k + 1L
    idx_within[j] <- k
  }
  s <- h2; s[s==""] <- NA
  for (j in seq_along(s)) {
    if (!is.na(grp[j]) && (is.na(s[j]) || !nzchar(s[j]))) {
      pos <- ((idx_within[j] - 1L) %% length(metrics_order)) + 1L
      s[j] <- metrics_order[pos]
    }
  }
  
  # Final column names: globals keep h2; allele-block columns become "ALLELE__METRIC"
  globals <- c("Pos","Peptide","ID","Target","Ave","NB")
  make_name <- function(g, ss) if (is.na(g) || ss %in% globals) ss else paste0(g, "__", ss)
  coln <- mapply(make_name, grp, s, USE.NAMES = FALSE)
  coln <- make.unique(coln, sep = "_dup")
  
  dat <- raw[-c(1,2), , drop = FALSE]
  names(dat) <- coln
  
  # Pivot only allele columns; globals are left as-is
  allele_cols <- grep("__", names(dat), value = TRUE)
  long <- dat %>%
    pivot_longer(cols = all_of(allele_cols),
                 names_to = c("Allele","Metric"),
                 names_sep = "__",
                 values_to = "value") %>%
    pivot_wider(names_from = Metric, values_from = value)
  
  # Coerce types
  suppressWarnings({
    long <- long %>%
      mutate(
        Pos      = parse_number(Pos),
        Inverted = parse_number(Inverted),
        Score    = parse_number(Score),
        Rank     = parse_number(Rank),
        Score_BA = parse_number(Score_BA),
        `nM`     = parse_number(`nM`),
        Rank_BA  = parse_number(Rank_BA),
        Ave      = parse_number(Ave),
        NB       = parse_number(NB)
      )
  })
  
  # Binder call from EL %Rank
  long %>%
    mutate(
      Binder = case_when(
        !is.na(Rank) & Rank <= strong_thr ~ "SB",
        !is.na(Rank) & Rank <= weak_thr   ~ "WB",
        TRUE                              ~ "NB"
      )
    ) %>%
    arrange(Allele, Rank, desc(Score))
}


write_fasta_with_genes <- function(sequences, genes, file) {
  # ensure equal length
  stopifnot(length(sequences) == length(genes))
  
  # make gene names unique if duplicated
  unique_genes <- make.unique(genes, sep = "_")
  
  con <- file(file, "w")
  for (i in seq_along(sequences)) {
    cat(paste0(">", unique_genes[i], "\n"), file = con)
    cat(sequences[i], "\n", file = con)
  }
  close(con)
}
