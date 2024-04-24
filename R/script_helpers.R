# Functions to facillitate Proteomics and Phosphoproteomics script execution

# Fix incorrect combiPeptData gene names ----
#' Fix incorrect gene names in Pescal Excel file
#'
#' Fixes incorrect data of pre-processed mass spectrometry Excel file:
#' combiPeptData columns 25 and 29.
#'
#' @param df.combi \code{\link{data.frame}}. combiPeptData.
#' @param organism String. One of 'human', 'mouse', 'rat', or 'pig'. Selects
#' database of Uniprot names.
#' @param pescal.xlsx String. The path ("directory/name") of the Pescal Excel
#' output.
#' @param fixed.xlsx String. Name for the fixed output Excel file.
#'
#' @return \code{\link{data.frame}}. The corrected combiPeptData. Additionally,
#' saves a copy of the Pescal Excel file with the corrected combiPeptData sheet
#' (with the name of `fixed.xlsx`).
#'
#' @export
fix_combiPeptData <- function(
    df.combi = df.combi,
    organism = organism_dbs,
    pescal.xlsx = pescal.output.file,
    fixed.xlsx = df.combi.fixed.xlsx
) {
  require(tidyverse)

  # Name columns
  colnames(df.combi)[30] <- "uniprot"
  colnames(df.combi)[29] <- "genes"
  colnames(df.combi)[25] <- "genes_mod"

  # Get uniprot data
  if (organism == "human") {
    df.uniprot <- protools2::uniprot_human
  } else if (organism == "mouse") {
    df.uniprot <- protools2::uniprot_mouse
  } else if (organism == "rat") {
    df.uniprot <- protools2::uniprot_rat
  } else if (organism == "pig") {
    df.uniprot <- protools2::uniprot_pig
  } else {
    stop("Error: `organism` must be one of 'human', 'mouse', 'rat', or 'pig'")
  }

  # Fix genes column
  df1 <- df.combi %>%
    dplyr::mutate(uniprot_split = stringr::str_split(uniprot, "; ")) %>%
    tidyr::unnest(uniprot_split) %>%
    dplyr::mutate(uniprot_split = stringr::str_replace_all(uniprot_split, ";", ""))

  df2 <- df.uniprot %>%
    dplyr::filter(Entry %in% df1$uniprot_split) %>%
    dplyr::select(Entry, Gene.Names..primary.)

  df3 <- df1 %>%
    dplyr::left_join(df2, by = c("uniprot_split" = "Entry")) %>%
    dplyr::group_by(db_id) %>%
    dplyr::summarise(genes_new = paste(Gene.Names..primary., collapse = "; ")) %>%
    dplyr::mutate(genes_new = paste0(genes_new, ";"))

  df4 <- df.combi %>%
    dplyr::left_join(df3, by = "db_id")

  # Fix genes_mod column
  df5 <- df4 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      genes_mod_vals = list(stringr::str_split_1(genes_mod, ";")),
      genes_mod_keep = list(unlist(stringr::str_extract_all(genes_mod_vals, "\\(.*?\\)"))),
      genes_new_split = list(stringr::str_replace_all(stringr::str_split_1(genes_new, "; "), ";", "")),
      genes_mod_reinsert = list(paste0(genes_new_split, genes_mod_keep)),
      genes_mod_recomb = paste0(genes_mod_reinsert, collapse = ";"),
      genes_mod_recomb = paste0(genes_mod_recomb, ";")
    ) %>%
    dplyr::mutate(
      genes = genes_new,
      genes_mod = genes_mod_recomb,
      genes_new = NULL,
      genes_mod_vals = NULL,
      genes_mod_keep = NULL,
      genes_new_split = NULL,
      genes_mod_reinsert = NULL,
      genes_mod_recomb = NULL
    ) %>%
    dplyr::rename(
      ...25 = genes_mod,
      ...29 = genes,
      ...30 = uniprot
    )

  # Fix incorrect entries in Excel file
  # Load existing Pescal Excel
  wb <- openxlsx::loadWorkbook(pescal.xlsx)
  # Overwrite combiPeptData sheet
  openxlsx::writeData(wb, sheet = "combiPeptData", x = df5, colNames = TRUE, rowNames = FALSE)
  # Save changes to new Excel file
  openxlsx::saveWorkbook(wb, fixed.xlsx, overwrite = TRUE)

  rm(df1, df2, df3, df4, wb)
  return(df5)
}


# Merge function ----
#' Merge multiple `data.table`s
#'
#' Merges multiple \code{\link{data.table}}s.
#'
#' @param dt_list \code{\link{list}} of \code{\link{data.table}}s.
#' @param by String/vector of strings. The column(s) to merge on.
#' @param all Logical. Default `TRUE` is an outer join; `FALSE` is an inner join.
#' @param sort Logical. To sort on the `by` columns. Default `FALSE`.
#'
#' @return A \code{\link{data.table}}. The result of the merge.
#' @export
#'
#' @examples
#' DT <- data.table(x=rep(c("b","a","c"),each=3), y=c(1,3,6), v=1:9)
#' A <- data.table(x=c("c","b"), v=8:7, foo=c(4,2))
#'
#' multi_DT <- mergeDTs(
#'   list(DT, A),
#'   by = "x",
#'   all = TRUE  # Outer join
#' )
mergeDTs <- function(dt_list, by = NULL, all = TRUE, sort = FALSE) {
  Reduce(
    function(...) {
      merge(..., by = by, all = all, sort = sort)
    },
    dt_list
  )
}


# Edited `protools2::normalize_areas_return_ppindex()` ----
# Phosphoproteomics.Rmd:
# Original function mistakenly used object name not of parameter, but of object in
# global environment. This meant that `df.combi` was using original from global scope not from argument
normalize_areas_return_ppindex_edit <- function(
    pescal_output_file,
    delta_score_cut_off = 0
) {
  # Set delta_score_cut_off to low (say 1) for proteomics data,
  # High (say 15) for phosphoproteomics data
  library(foreach)
  library(doParallel)
  suppressMessages(
    df.areas <- readxl::read_excel(pescal_output_file, "output_areas")
  )
  colnames(df.areas) <- gsub("-", ".", colnames(df.areas), fixed = T)
  suppressMessages(
    df.combi <- readxl::read_excel(pescal_output_file, "combiPeptData")  # Original `pescal.output.file` object is accessed from global scope
  )

  # select peptides above the delta_score_cut_off
  df.combi <- subset(df.combi, df.combi$max_delta_score > delta_score_cut_off)
  peptides <- unique(unlist(df.combi[, 25]))  # 'sites'
  df.areas <- df.areas[df.areas$db_id %in% df.combi$db_id, ]
  cols <- colnames(dplyr::select_if(df.areas, is.numeric))
  df.areas.n <- data.frame(
    ids = df.areas$db_id,
    scale(df.areas[, cols], center = F, scale = colSums(df.areas[, cols]))
  )

  cores = detectCores()
  cl <- makeCluster(cores[1] - 1)  # not to overload your computer
  registerDoParallel(cl)
  t1 <- Sys.time()
  df <- foreach(p = peptides, .combine = "rbind") %dopar%
    {
      ids <- na.omit(df.combi[df.combi[, 25] == p, ]$db_id)
      apply(df.areas.n[df.areas.n$ids %in% ids, cols], 2, sum)
    }
  stopCluster(cl)
  df.norm <- data.frame(sites = peptides, df * 1e+06)
  rownames(df.norm) <- df.norm$sites
  df.norm[df.norm == 0] <- NA

  df.norm.log2.centered <- data.frame(
    sites = peptides,
    scale(log2(df.norm[, cols]), scale = F)
  )

  df.norm.log2.centered.scaled <- data.frame(
    sites = peptides,
    scale(log2(df.norm[, cols]))
  )

  # Alternative na imputation
  # Centred & scaled
  df.norm.log2.centered.scaled.na.imputed <- df.norm.log2.centered.scaled
  df.norm.log2.centered.scaled.na.imputed1 <- df.norm.log2.centered.scaled
  df.norm.log2.centered.scaled.na.imputed2 <- df.norm.log2.centered.scaled

  df.norm.log2.centered.scaled.na.imputed[cols] <- lapply(
    df.norm.log2.centered.scaled.na.imputed[cols], function(x){
      replace(x, is.na(x), min(x, na.rm = TRUE) -1) # Correct NA imputation
    }
  )

  df.norm.log2.centered.scaled.na.imputed1[
    is.na(df.norm.log2.centered.scaled.na.imputed1)
  ] <- min(df.norm.log2.centered.scaled.na.imputed1[,cols], na.rm = T) / 5

  df.norm.log2.centered.scaled.na.imputed2[cols] <- lapply(
    df.norm.log2.centered.scaled.na.imputed2[cols], function(x){
      replace(x, is.na(x), min(x, na.rm = TRUE)) # Correct NA imputation
    }
  )

  # Scaled
  df.norm.log2.centered.na.imputed <- df.norm.log2.centered
  df.norm.log2.centered.na.imputed[cols] <- lapply(
    df.norm.log2.centered.na.imputed[cols], function(x){
      replace(x, is.na(x), min(x, na.rm = TRUE) -1) # Correct NA imputation
    }
  )

  # Previous na imputation
  # df.norm.log2.centered.scaled.na.imputed <- df.norm.log2.centered.scaled
  # df.norm.log2.centered.scaled.na.imputed[
  #   is.na(df.norm.log2.centered.scaled.na.imputed)
  # ] <- min(df.norm.log2.centered.scaled.na.imputed[, cols], na.rm = T) / 5

  rownames(df.norm.log2.centered) <- df.norm.log2.centered$sites
  rownames(df.norm.log2.centered.scaled) <- df.norm.log2.centered.scaled$sites

  return(list(
    normalized.data = df.norm,
    normalized.plus.log2.cent.data = df.norm.log2.centered,
    normalized.plus.log2.cent.scaled.data = df.norm.log2.centered.scaled,
    df.norm.log2.centered.scaled.na.imputed = df.norm.log2.centered.scaled.na.imputed,
    df.norm.log2.centered.na.imputed=df.norm.log2.centered.na.imputed,
    df.norm.log2.centered.scaled.na.imputed1=df.norm.log2.centered.scaled.na.imputed1,
    df.norm.log2.centered.scaled.na.imputed2=df.norm.log2.centered.scaled.na.imputed2
  ))
}


# Edited `protools2::normalize_areas_return_protein_groups()` ----
# Proteomics.Rmd:
# Original function mistakenly used object name not of parameter, but of object in
# global environment. This meant that `df.combi` was using original from global scope not from argument
normalize_areas_return_protein_groups_edit <- function(
    pescal_output_file,
    mascot.score.cut.off = 50,
    n.peptide.cut.off = 1
) {
  suppressMessages(
    df.areas <- data.frame(readxl::read_excel(pescal_output_file, "output_areas"))
  )
  colnames(df.areas) <- gsub("-", ".", colnames(df.areas), fixed = T)
  suppressMessages(
    df.combi <- data.frame(readxl::read_excel(pescal_output_file, "combiPeptData"))  # Original `pescal.output.file` object is accessed from global scope
  )

  # normalise areas
  cols <- colnames(dplyr::select_if(df.areas, is.numeric))
  df.areas.n <- data.frame(
    ids = df.areas$db_id,
    scale(df.areas[, cols], center = F, scale = colSums(df.areas[, cols]))
  )

  # find protein groups
  protein.groups <- na.omit(unique(unlist(df.combi[, 29])))  # 29 = genes
  n.protein.groups <- length(protein.groups)

  # group peptides by protein group
  cores = detectCores()
  cl <- makeCluster(cores[1] - 1)  # not to overload your computer
  registerDoParallel(cl)
  t1 <- Sys.time()
  df <- foreach(p = protein.groups, .combine = "rbind") %dopar%
    {
      dfx <- df.combi[df.combi[, 29] == p, ]
      ids <- na.omit(dfx$db_id)
      best.mascot.score <- max(dfx$max_scr, na.rm = T)
      protein.name <- dfx$protein[1]
      acc <- na.omit(dfx$acc_no)[1]
      uniprot.id <- na.omit(dfx[1, 30])[1]
      n.peptides <- length(ids)
      nPSMs <- na.omit(dfx[, "N_peptides"])
      c(
        protein.group = p,
        apply(df.areas.n[df.areas.n$ids %in% ids, cols], 2, sum),
        best.mascot.score = best.mascot.score,
        n.peptides = n.peptides,
        n.psm = sum(nPSMs),
        acc = acc,
        uniprot.id = uniprot.id,
        protein.name = protein.name
      )
    }
  stopCluster(cl)

  write.csv(df, "temp.csv")
  x <- read.csv("temp.csv")
  x[x == 0] <- NA
  df.norm <- data.frame(protein.group = x$protein.group, x[, cols] * 1e+06)

  df.norm.log2.centered <- data.frame(
    protein.group = protein.groups,
    scale(log2(df.norm[, cols]), scale = F)
  )

  df.norm.log2.centered.scaled <- data.frame(
    protein.group = protein.groups,
    scale(log2(df.norm[, cols]))
  )

  df.norm.log2.centered.scaled.na.imputed <- df.norm.log2.centered.scaled

  # df.norm.log2.centered.scaled.na.imputed[  # Previous na imputation
  #   is.na(df.norm.log2.centered.scaled.na.imputed)
  # ] <- min(df.norm.log2.centered.scaled.na.imputed[, cols], na.rm = T) / 5

  df.norm.log2.centered.scaled.na.imputed[cols] <- lapply(  # New na imputation
    df.norm.log2.centered.scaled.na.imputed[cols], function(x){
      replace(x, is.na(x), min(x, na.rm = TRUE) -1) # Correct NA imputation
    }
  )

  rownames(df.norm.log2.centered) <- df.norm.log2.centered$protein.group
  rownames(df.norm.log2.centered.scaled) <- df.norm.log2.centered.scaled$protein.group
  rownames(df.norm) <- df.norm.log2.centered$protein.group
  rownames(df.norm.log2.centered.scaled.na.imputed) <- df.norm.log2.centered.scaled.na.imputed$protein.group

  xx <- x[
    x$best.mascot.score > mascot.score.cut.off &
      x$n.peptides > n.peptide.cut.off, ]

  selected.prot.groups <- xx$protein.group

  cc <- c(
    "protein.group", "best.mascot.score", "n.peptides",
    "n.psm", "acc", "uniprot.id", "protein.name"
  )

  df.norm <- merge.data.frame(df.norm, x[, cc], by = "protein.group")
  df.norm.log2.centered <- merge.data.frame(
    df.norm.log2.centered,
    x[, cc],
    by = "protein.group"
  )
  df.norm.log2.centered.scaled <- merge.data.frame(
    df.norm.log2.centered.scaled,
    x[, cc],
    by = "protein.group"
  )
  df.norm.log2.centered.scaled.na.imputed <- merge.data.frame(
    df.norm.log2.centered.scaled.na.imputed,
    x[, cc],
    by = "protein.group"
  )

  return(list(
    normalized.data = df.norm[df.norm$protein.group %in% selected.prot.groups, ],
    normalized.plus.log2.cent.data = df.norm.log2.centered[
      df.norm.log2.centered$protein.group %in% selected.prot.groups, ],
    normalized.plus.log2.cent.scaled.data = df.norm.log2.centered.scaled[
      df.norm.log2.centered.scaled$protein.group %in% selected.prot.groups, ],
    df.norm.log2.centered.scaled.na.imputed = df.norm.log2.centered.scaled.na.imputed[
      df.norm.log2.centered.scaled.na.imputed$protein.group %in% selected.prot.groups, ]
  ))
}


# Edited `protools2::remove_outlier_samples_from_dataset()` ----
# `protools2::remove_outlier_samples_from_dataset()` returns a NULL object
# if no outliers are present, edited function instead returns the original
# dataframe.

#' Remove outliers from a dataframe.
#'
#' Takes a dataframe in wide format, calculates medians for each numeric
#'   column, and uses `boxplot()$out` values to identify outliers for removal.
#'
#' @param df A `data.frame`.
#' @param plot Boolean. If TRUE (default) then the boxplot is produced. If
#'   FALSE then no plot is produced.
#'
#' @return Returns the given dataset with columns identified as outliers
#'   removed. If no outliers are present in the dataset, then the original
#'   dataset is returned.
#'
#' @seealso \code{\link{boxplot}}
#'
#' @export
remove_outlier_samples_from_dataset_edit <- function(df, plot=TRUE) {
  cols <- colnames(dplyr::select_if(df, is.numeric))
  cms <- apply(df[, cols], 2, median)
  if (plot) {
    OutVals <- boxplot(cms)$out  # main = "Outliers"
  } else {
    OutVals <- boxplot(cms, plot=FALSE)$out  # main = "Outliers"
  }
  samples.to.be.removed <- names(OutVals)
  if (is.null(samples.to.be.removed)) {
    return(df)
  } else {
    df_res <- df[, -which(names(df) %in% samples.to.be.removed)]
    return(df_res)
  }
}


# Edited `protools2::summary.qual.data()` ----
# To handle non-equal replicates in the samples
summary.qual.data_edit <- function(df.combi, df.design) {
  df.mods <- data.frame(table(df.combi$pep_mod))
  df.mods.phospho.st <- df.mods[grepl("ST", df.mods$Var1), ]
  sites <- df.combi[, 25]
  phospho.s.unique <- format(
    nrow(sites[grepl("S", sites[[1]], fixed = T), ]), big.mark = ","
  )
  phospho.t.unique <- format(
    nrow(sites[grepl("(T", sites[[1]], fixed = T), ]), big.mark = ","
  )
  phospho.y.unique <- format(
    nrow(sites[grepl("(Y", sites[[1]], fixed = T), ]), big.mark = ","
  )
  n.st.peptides <- sum(df.mods.phospho.st$Freq)
  df.mods.phospho.y <- df.mods[grepl("(Y)", df.mods$Var1, fixed = T), ]
  n.y.peptides <- format(sum(df.mods.phospho.y$Freq), big.mark = ",")
  n.peptides.unique <- format(length(unique(unlist(sites))), big.mark = ",")
  n.peptides.total <- format(nrow(df.combi), big.mark = ",")
  nconditions <- length(unique(df.design$condition))
  nruns <- nrow(df.design)
  df.phospho <- subset(df.combi, grepl("Phospho", df.combi$pep_mod, fixed = T))
  total.phosphoproteins <- format(length(unique(df.phospho$acc_no)), big.mark = ",")
  total.proteins <- format(length(unique(df.combi$acc_no)), big.mark = ",")
  summary.text <- paste(
    "The experiment identified", n.peptides.total,
    "peptides, of which", n.peptides.unique, "were unique. There were",
    phospho.s.unique, "unique phosphoserine,", phospho.t.unique,
    "phosphothreonine and", phospho.y.unique, "phosphotyrosine sites identified and quantified in the experiment.",
    "A total of", total.proteins, "proteins were identified in the experiment of which",
    total.phosphoproteins, "were phosphorylated.", " \n",
    "The experiment compared", nconditions, "conditions in",
    if ((nruns / nconditions) %% 1 == 0) nruns / nconditions else paste("a total of", nruns),
    "replicates, requiring", nruns, "LC-MS/MS runs."
  )
}


# Edited `protools2::summary.qual.data.proteomics()` ----
# Edited as `nruns` was defined after use
summary.qual.data.proteomics_edit <- function(
    df.norm,
    df.design,
    replicates = 0,
    conditions = 0
) {
  protein.groups <- unique(unlist(df.norm$protein.group))
  n.protein.groups <- format(length(protein.groups), big.mark = ",")
  n.peptides.unique <- format(sum(df.norm$n.peptides), big.mark = ",")
  n.peptides.total <- format(sum(df.norm$n.psm), big.mark = ",")
  total.proteins <- format(length(unique(df.norm$acc)), big.mark = ",")
  nruns <- replicates * conditions
  if (replicates == 0) {
    conditions <- length(unique(df.design$condition))
    replicates <- nruns/conditions
  }
  data.points <- format(nruns * length(protein.groups), big.mark = ",")
  summary.text <- paste(
    "The experiment produced", n.peptides.total,
    "peptides spectral matches, of which", n.peptides.unique,
    "were unique peptides, and belonging to", n.protein.groups,
    "protein groups.", " \n", "The experiment compared",
    conditions, "conditions in", replicates, "replicates, requiring",
    nruns, "LC-MS/MS runs and producing", data.points, "quantitative datapoints."
  )
}


# Edited `protools2::pca.plot()` ----
pca.plot_edit <- function(df, df.design, colorfactor = "", shapefactor = "", legend_position = NULL, title = "") {
  df <- dplyr::select_if(df, is.numeric)
  res.pca <- prcomp(df, scale = F)
  df.pca <- data.frame(res.pca$rotation)
  df.pca <- data.frame(heading = rownames(df.pca), df.pca)
  df.pca <- merge.data.frame(df.pca, df.design, by = "heading")
  ss <- summary(res.pca)
  df.sum <- data.frame(ss$importance)
  pc1 <- df.sum$PC1[2]
  pc2 <- df.sum$PC2[2]
  pc3 <- df.sum$PC3[2]
  # f <- round(nrow(df.pca)/20)
  f <- nrow(df.pca)  # Above results in error as too few to create scales.
  plot1 <- ggplot(
    df.pca,
    aes_string(x = "PC1", y = "PC2", color = colorfactor, shape = shapefactor)
  ) +
    geom_point() +
    scale_color_manual(values = rep(mycolors()$C23, f)) +
    scale_shape_manual(values = rep(c(1:20), f)) +
    theme_bw() +
    labs(
      title = title,
      x = paste("PC1,", round(pc1, digits = 2) * 100, "%"),
      y = paste("PC2,", round(pc2, digits = 2) * 100, "%")
    )
  if (!is.null(legend_position)) {
    plot1 <- plot1 + theme(
      legend.position = legend_position,
      legend.text = element_text(
        size = 8,
        lineheight = .8
      )
    )
    legend1 <- cowplot::get_legend(plot1)
    plot1 <- plot1 + theme(legend.position = "none")
  }
  plot2 <- ggplot(
    df.pca,
    aes_string(x = "PC2", y = "PC3", color = colorfactor, shape = shapefactor)
  ) +
    geom_point() +
    scale_color_manual(values = rep(mycolors()$C23, f)) +
    scale_shape_manual(values = rep(c(1:20), f)) + theme_bw() +
    labs(
      title = "",
      x = paste("PC2,", round(pc2, digits = 2) * 100, "%"),
      y = paste("PC3,", round(pc3, digits = 2) * 100, "%")
    )
  if (!is.null(legend_position)) {
    plot2 <- plot2 + theme(
      legend.position = legend_position,
      legend.text = element_text(
        size = 8,
        lineheight = .8
      )
    )
    legend2 <- cowplot::get_legend(plot2)
    plot2 <- plot2 + theme(legend.position = "none")
  }
  return(
    if (is.null(legend_position)) {
      cowplot::plot_grid(plot1, plot2)
    } else if (legend_position == "bottom" | legend_position == "top") {
      cowplot::plot_grid(plot1, plot2, legend1, nrow = 2, ncol = 2)
    } else if (legend_position == "left" | legend_position == "right") {
      cowplot::plot_grid(plot1, plot2, legend1, nrow = 1, ncol = 3)
    }
  )
}


# Edited `protools2::identify_differences_in_comparison_plus_volcano()` ----
identify_differences_in_comparison_plus_volcano_edit <- function (df, fold.cutoff = 0.5, qval.cutoff = 0.05, graph.header = "")
{
  colnames(df) <- c("sites", "fold", "pvalue", "qvalue")
  df.up <- subset(df, df$qvalue < qval.cutoff & df$fold > fold.cutoff)
  df.do <- subset(df, df$qvalue < qval.cutoff & df$fold < (-fold.cutoff))
  df.up <- df.up[order(-df.up$fold), ]
  df.do <- df.do[order(df.do$fold), ]
  suppressMessages(
    plot.pval.dist <- ggplot(df, aes(x = pvalue)) +
      geom_histogram(
        # fill = "white",
        alpha = 0.7,
        color = "black"
      ) +
      labs(title = "Distrubtion of pvalues") +
      theme_bw()
  )
  suppressMessages(
    plot.qval.dist <- ggplot(df, aes(x = qvalue)) +
      geom_histogram(
        # fill = "white",
        alpha = 0.7,
        color = "black"
      ) +
      labs(title = "Distrubtion of qvalues (BH adjusted pvalues)") +
      theme_bw()
  )
  plot.volcano.1 <- ggplot(df, aes(x = fold, y = -log10(qvalue))) +
    geom_point(size = 1) +
    geom_point(data = df.up, aes(x = fold, y = -log10(qvalue)), color = "red") +
    geom_point(data = df.do, aes(x = fold, y = -log10(qvalue)), color = "blue") +
    labs(
      title = graph.header,
      subtitle = paste0("Increased = ", nrow(df.up), "; Decreased = ", nrow(df.do))
    )
  return(
    list(
      df.increased = df.up,
      df.decreased = df.do,
      volcanoplot = plot.volcano.1,
      pvalue.distributions = cowplot::plot_grid(plot.pval.dist, plot.qval.dist)
    )
  )
}


# Edited `protools2::identify_differences_by_pvalue_in_comparison_plus_volcano()` ----
# Original incorrectly displays duplicated pvalue distribution, so removed repeated plot
identify_differences_by_pvalue_in_comparison_plus_volcano_edit <- function(df, fold.cutoff = 0.5, pval.cutoff = 0.05, graph.header = "") {
  colnames(df) <- c("sites", "fold", "pvalue", "qvalue")
  df.up <- subset(df, df$pvalue < pval.cutoff & df$fold >
                    fold.cutoff)
  df.do <- subset(df, df$pvalue < pval.cutoff & df$fold <
                    (-fold.cutoff))
  df.up <- df.up[order(-df.up$fold), ]
  df.do <- df.do[order(df.do$fold), ]
  suppressMessages(plot.pval.dist <- ggplot(df, aes(x = pvalue)) +
                     geom_histogram(
                       # fill = "white",
                       color = "black",
                       alpha = 0.7
                     ) +
                     labs(title = "Distrubtion of pvalues") +
                     theme_bw())
  suppressMessages(plot.qval.dist <- ggplot(df, aes(x = qvalue)) +
                     geom_histogram(
                       # fill = "white",
                       color = "black",
                       alpha = 0.7
                     ) +
                     labs(title = "Distrubtion of pvalues (BH adjusted pvalues)") +
                     theme_bw())
  plot.volcano.1 <- ggplot(df, aes(x = fold, y = -log10(pvalue))) +
    geom_point(size = 1) +
    geom_point(data = df.up, aes(
      x = fold,
      y = -log10(pvalue)
    ), color = "red") +
    geom_point(
      data = df.do,
      aes(x = fold, y = -log10(pvalue)), color = "blue"
    ) +
    labs(title = graph.header, subtitle = paste0(
      "Increased = ",
      nrow(df.up), "; Decreased = ", nrow(df.do)
    ))
  return(
    list(
      df.increased = df.up, df.decreased = df.do,
      volcanoplot = plot.volcano.1,
      pvalue.distributions = cowplot::plot_grid(plot.pval.dist, ggplot() + theme_void())
      # pvalue.distributions = cowplot::plot_grid(plot.pval.dist, plot.qval.dist)
    )
  )
}


# Edited `protools2::barplot.top.peptides()` ----
# Edited to change plot titles
barplot.top.peptides_edit <- function(
    df.up,
    df.do,
    graph.heading = "",
    context = "peptides",
    subtitle = ""
) {
  df.up$effect <- "Increased"
  df.do$effect <- "Decreased"
  df.up.down <- na.omit(
    rbind.data.frame(df.up[1:25, ], df.do[1:25, ])
  )
  # a <- df.up.down$sites  # The following 3 lines of code are replaced with `stringr` method below
  # a <- ifelse(nchar(a) > 20, paste0(strtrim(a, 20), "..."), a)
  # df.up.down$sites <- a
  df.up.down$sites <- stringr::str_trunc(df.up.down$sites, width = 30, side = "right")
  barplot.top.up <- ggplot(
    df.up.down,
    aes(x = reorder(sites, fold), y = fold)
  ) +
    geom_bar(stat = "identity", aes(fill = effect)) +
    scale_fill_manual(values = c("skyblue", "red")) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    labs(
      title = paste0("Top increased / decreased ", context," in ", graph.heading),
      y = "log2(fold)",
      x = "",
      subtitle = subtitle
    )
}


# Edited `protools2::plot.kinase.relationships()` ----
# `reshape2` is deprecated
plot.kinase.relationships_edit <- function(ksea.data, graph.title, pval.cutoff) {
  library(ggrepel)
  library(igraph)

  kinases <- ksea.data[ksea.data$pvalues < pval.cutoff, ]$kinase_dataset
  if (length(kinases) == 0) {
    return(0)
  }
  set.seed(123)
  pth1 <- character()
  pth2 <- character()
  pval1 <- numeric()
  pval2 <- numeric()
  zscores <- numeric()
  w <- numeric()
  g <- character()
  counts <- numeric()
  i <- 1
  for (k1 in kinases) {
    count <- ksea.data[ksea.data$kinase_dataset == k1, "m"]
    pvalue <- ksea.data[ksea.data$kinase_dataset == k1, "pvalues"]
    zscore <- ksea.data[ksea.data$kinase_dataset == k1, "zscores"]
    genes.1 <- unlist(strsplit(
      ksea.data[ksea.data$kinase_dataset == k1, "sites"], ";"
    ))
    for (k2 in kinases) {
      genes.2 <- unlist(
        strsplit(ksea.data[ksea.data$kinase_dataset == k2, "sites"], ";")
      )
      common <- intersect(genes.1, genes.2)
      if (length(common) > 0) {
        pth1[i] <- k1
        pth2[i] <- k2
        g[i] <- paste0(common, collapse = ";")
        w[i] <- length(common)
        pval1[i] <- ksea.data[ksea.data$kinase_dataset == k2, "pvalues"]
        pval2[i] <- ksea.data[ksea.data$kinase_dataset == k1, "pvalues"]
        zscores[i] <- ksea.data[ksea.data$kinase_dataset == k1, "zscores"]
        counts[i] <- ksea.data[ksea.data$kinase_dataset == k1, "m"]
        i <- i + 1
      }
    }
  }
  df.network <- data.frame(
    node.1 = pth1, node.2 = pth2, weight = w,
    pvalue.1 = pval1, pvalue.2 = pval2, zscores = zscores,
    counts = counts
  ) # POSSIBLY CHECK `if (nrow(df.network) < 2) { return(0) }`
  df.network <- df.network[order(df.network$pvalue.1), ]
  df.s <- subset(df.network, df.network$weight > 4) # NOTE: weight subset can produce no results
  if (nrow(df.s) == 0) {
    return(0)
  }
  xx <- reshape2::dcast(df.s, node.1 ~ node.2, value.var = "weight")
  xx[is.na(xx)] <- 0
  rownames(xx) <- xx$node.1
  caught.inc <- graph.incidence(xx[, 2:ncol(xx)], weighted = TRUE)
  obs.parties.all <- bipartite.projection(caught.inc)[[1]]
  obs.spp.all <- bipartite.projection(caught.inc)[[2]]
  fr.all <- layout_with_graphopt(obs.spp.all)
  fr.all.df <- as.data.frame(fr.all)
  fr.all.df$species <- colnames(xx[, 2:ncol(xx)]) # RESULT WAS `NULL` (Now works ? after `ksea_edit()`)
  if (is.null(colnames(xx[, 2:ncol(xx)]))) {  # possibly `NULL` if nrow() < 2 ?
    fr.all.df$species <- colnames(xx[2:ncol(xx)])  # As above now works, this is no longer needed (still is needed)
  }
  # EMPTY DATAFRAME `g` (0 rows, columns = "from", "to"):
  g <- get.data.frame(obs.spp.all, what = "edges") # WAS EMPTY DATAFRAME (Now works ? after `ksea_edit()`)
  # Attempt to recreate result of above line:
  if (nrow(g) == 0) {
    g <- df.s[, 1:3]  # Added - NO LONGER NEEDED (still is needed)
    colnames(g) <- c("from", "to", "weight")  # Added - NO LONGER NEEDED (still is needed)
  }
  g$from.x <- fr.all.df$V1[match(g$from, fr.all.df$species)]
  g$from.y <- fr.all.df$V2[match(g$from, fr.all.df$species)]
  g$to.x <- fr.all.df$V1[match(g$to, fr.all.df$species)]
  g$to.y <- fr.all.df$V2[match(g$to, fr.all.df$species)]
  g$pvalue.1 <- df.network$pvalue.1[match(g$from, df.network$node.1)]
  g$pvalue.2 <- df.network$pvalue.2[match(g$from, df.network$node.2)]
  g$weight <- df.network$weight[match(g$from, df.network$node.1)]
  fr.all.df$pvalue.1 <- df.network$pvalue.1[match(
    fr.all.df$species, df.network$node.1
  )]
  fr.all.df$pvalue.2 <- df.network$pvalue.2[match(
    fr.all.df$species, df.network$node.2
  )]
  fr.all.df$E <- df.network$zscores[match(
    fr.all.df$species, df.network$node.1 #
  )]
  fr.all.df$counts <- df.network$counts[match(
    fr.all.df$species, df.network$node.2
  )]
  pplot <- ggplot() +
    geom_segment(
      data = g, aes(
        x = from.x, xend = to.x, y = from.y, yend = to.y, size = weight / 10
      ),
      colour = "black", linetype = 2
    ) +
    geom_point(
      data = fr.all.df,
      aes(x = V1, y = V2, fill = E, size = counts), shape = 21
    ) +
    geom_label_repel(
      data = fr.all.df,
      aes(x = V1, y = V2, label = species)
    ) +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", midpoint = 0
    ) +
    scale_x_continuous(expand = expansion(add = 10)) + # c(0, 1)) +  # `expansion()` to provide padding
    scale_y_continuous(expand = expansion(add = 10)) + # c(0, 1)) +
    labs(y = "") + # To prevent overlapping of labels and title
    theme_bw() +
    theme(
      axis.text.x = element_blank(), axis.text.y = element_blank(),
      axis.ticks = element_blank(), axis.title.x = element_blank(),
      axis.title.y = element_text(size = 30), panel.background = element_blank(), # axis.title.y modified
      panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), plot.background = element_blank()
    ) +
    ggtitle(graph.title) +
    # labs(caption = stringr::str_wrap(
    #   "Edge weights are proportional to common downstream targets between kinases. Node colors are proportional to z-score of enrichment.",
    #   width = 75
    # ))
    labs(
      caption = paste(
        "Edges indicate common neighbours.",
        "Edge weights are proportional to common downstream targets between kinases.",
        "Node colors are proportional to z-score of enrichment.",
        sep = "\n"
      )
    )
  pplot
  return(pplot)
}


# Edited `protools2::pathway_enrichment()` ----
# Edited to adjust colour scales
pathway_enrichment_edit <- function (
    increased.peptides,
    decreased.peptides,
    background.data,
    prot_dbs = c("kegg", "hallmark.genes", "nci", "process"),
    graph.heading = "",
    is.ksea = FALSE
) {
  library(foreach)
  library(doParallel)
  background.list <- background.data$protein
  cores = detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  t1 <- Sys.time()
  enrich.up <- foreach(db = prot_dbs, .combine = "rbind") %dopar%
    {
      e <- protools2::enrichment.from.list(
        list.of.peptides = increased.peptides,
        background.list,
        prot_db = db,
        is.ksea = is.ksea
      )
      if (nrow(e) > 0) {
        e$prot_db <- db
        if (is.ksea) {
          e$pathway <- paste0(e$pathway, "_", db)
        }
        e
      }
    }
  stopCluster(cl)
  t2 <- Sys.time()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  t1 <- Sys.time()
  enrich.do <- foreach(db = prot_dbs, .combine = "rbind") %dopar%
    {
      e <- protools2::enrichment.from.list(
        list.of.peptides = decreased.peptides,
        background.list, prot_db = db, is.ksea = is.ksea
      )
      if (nrow(e) > 0) {
        e$prot_db <- db
        if (is.ksea) {
          e$pathway <- paste0(e$pathway, "_", db)
        }
        e
      }
    }
  stopCluster(cl)
  t2 <- Sys.time()
  if (length(enrich.do) > 0 & length(enrich.up) > 0) {
    enrich.up <- enrich.up[order(enrich.up$pvalues), ]
    a <- enrich.up$pathway
    a <- ifelse(nchar(a) > 50, paste0(strtrim(a, 50), "..."), a)
    enrich.up$pathway <- a
    enrich.do <- enrich.do[order(enrich.do$pvalues), ]
    a <- enrich.do$pathway
    a <- ifelse(nchar(a) > 50, paste0(strtrim(a, 50), "..."), a)
    enrich.do$pathway <- a
    enrich.do$effect <- "Decreased"
    enrich.up$effect <- "Increased"
    enrich.up.do <- rbind.data.frame(
      enrich.do[1:25, ], enrich.up[1:25, ]
    )
    diff.enrich <- merge.data.frame(
      enrich.do, enrich.up, by = "pathway", all = T
    )
    diff.enrich$enrichment.x[is.na(diff.enrich$enrichment.x)] <- 0
    diff.enrich$enrichment.y[is.na(diff.enrich$enrichment.y)] <- 0
    diff.enrich$pvalues.x[is.na(diff.enrich$pvalues.x)] <- 1
    diff.enrich$pvalues.y[is.na(diff.enrich$pvalues.y)] <- 1
    diff.enrich$counts.x[is.na(diff.enrich$counts.x)] <- 0
    diff.enrich$counts.y[is.na(diff.enrich$counts.y)] <- 0
    diff.enrich$delta.enrichment.up.vs.down <- (
      diff.enrich$enrichment.y - diff.enrich$enrichment.x
    )
    diff.enrich$delta.pval.up.vs.down <- -log10(diff.enrich$pvalues.y) -
      (-log10(diff.enrich$pvalues.x))
    diff.enrich$delta.counts <- diff.enrich$counts.y - diff.enrich$counts.x
    diff.enrich <- diff.enrich[order(-diff.enrich$delta.enrichment.up.vs.down), ]
    dif.e.up <- diff.enrich[1:15, ]
    diff.enrich <- diff.enrich[order(diff.enrich$delta.enrichment.up.vs.down), ]
    dif.e.do <- diff.enrich[1:15, ]
    dothis <- 0
    if (dothis == 1) {
      cores = detectCores()
      cl <- makeCluster(cores[1] - 1)
      registerDoParallel(cl)
      t1 <- Sys.time()
      enrich.combined <- foreach(db = prot_dbs, .combine = "rbind") %dopar%
        {
          e <- protools2::enrichment.from.list(
            list.of.peptides = c(increased.peptides, decreased.peptides),
            background.list,
            prot_db = db,
            is.ksea = is.ksea
          )
          if (nrow(e) > 0) {
            e$prot_db <- db
            e
          }
        }
      stopCluster(cl)
      t2 <- Sys.time()
      enrich.combined <- enrich.combined[order(enrich.combined$pvalues), ]
      a <- enrich.combined$pathway
      a <- ifelse(nchar(a) > 50, paste0(strtrim(a, 50), "..."), a)
      enrich.combined$pathway <- a
      enrich.combined <- enrich.combined[order(enrich.combined$pvalues), ]
      plot.pathways.enrichment.by.counts.combined <- ggplot(
        na.omit(enrich.combined),
        aes(x = counts, y = reorder(pathway, counts))
      ) +
        geom_point(aes(size = log2(enrichment), color = -log10(FDR))) +
        labs(
          title = "Pathway Enrichement Analysis",
          subtitle = "Combined analysis using increased and decreased peptides",
          size = "log2(E)", color = "-log10(q)",
          x = "# Proteins", y = "") +
        scale_color_gradient(
          low = "orange", high = "red",
          oob = scales::squish_infinite
        ) +  # na.value = "red") +
        facet_wrap(~effect, scales = "free") +
        theme_bw()
    }
    # Used in report
    plot.diff <- ggplot(
      rbind.data.frame(dif.e.do, dif.e.up),
      aes(
        x = delta.enrichment.up.vs.down,
        y = reorder(pathway, delta.enrichment.up.vs.down)
      )
    ) +
      geom_point(aes(size = abs(delta.counts), color = delta.pval.up.vs.down)) +
      labs(
        title = "Pathway Enrichement Analysis",
        subtitle = paste0("Changed in ", graph.heading),
        y = "", size = "delta \ncounts", color = "delta \npvalue"
      ) +
      scale_color_gradient2(
        low = "blue", mid = "yellow", high = "red",
        oob = scales::squish_infinite
      ) +  # na.value = "red") +
      # scale_colour_viridis(option = "turbo") +
      geom_vline(xintercept = 0, linetype = 2) +
      theme_bw()
    # Used in report
    plot.pathways.enrichment.by.counts <- ggplot(
      na.omit(enrich.up.do),
      aes(x = counts, y = reorder(pathway, counts))
    ) +
      geom_point(aes(size = log2(enrichment), color = -log10(FDR))) +
      labs(
        title = "Pathway Enrichement Analysis",
        subtitle = paste0("Changed in ", graph.heading),
        size = "log2(E)", color = "-log10(q)",
        x = "# Proteins", y = ""
      ) +
      scale_color_gradient(
        low = "orange", high = "red",
        oob = scales::squish_infinite
      ) +  # na.value = "red") +
      # scale_colour_viridis(option = "plasma") +
      facet_wrap(~effect, scales = "free") +
      theme_bw()
    plot.pathways.enrichment.by.enrichment <- ggplot(
      na.omit(enrich.up.do),
      aes(x = log2(enrichment), y = reorder(pathway, enrichment))
    ) +
      geom_point(aes(size = counts, color = -log10(FDR))) +
      labs(
        title = "Pathway Enrichement Analysis",
        subtitle = paste0("Changed in ", graph.heading),
        color = "-log10(q)", y = ""
      ) +
      scale_color_gradient(
        low = "orange", high = "red",
        oob = scales::squish_infinite
      ) +  # na.value = "red") +
      facet_wrap(~effect, scales = "free") +
      theme_bw()
    uniprot.data <- protools2::uniprot.names
    enrich.data <- enrich.up
    .look.at.proteins.in.pathways <- function(enrich.data) {
      mypathways <- list()
      j <- 1
      mypathways <- foreach(j = 1:10, .combine = "rbind") %do%
        {
          ontology <- enrich.data[j, ]$pathway
          prots <- unlist(strsplit(enrich.data[j, "proteins"], ";"))
          prot <- prots[1]
          dfxx <- foreach(prot = prots, .combine = "rbind") %do%
            {
              vv <- background.data[grepl(prot, background.data$protein, fixed = T), ]
              vv <- vv[order(vv$pvalues), ]
              vv[1, ]
            }
          prots <- unique(unlist(strsplit(dfxx$protein, ";")))
          df.names <- foreach(
            ppp = unlist(dfxx$protein), .combine = "rbind") %do% {
              prot <- unlist(strsplit(ppp, ";"))
              fff <- foreach(pp = prot, .combine = "rbind") %do%
                {
                  p <- strsplit(pp, "(", fixed = T)[[1]][1]
                  nn <- uniprot.data[grepl(p, uniprot.data$gene.name, fixed = T), ]
                  accs <- paste0(nn$acc, collapse = ";")
                  names <- paste0(nn$protein, collapse = ";")
                  data.frame(protein = ppp, accs, names)
                }
              data.frame(
                protein = ppp,
                acc = paste0(fff$accs, collapse = ";"),
                name = paste0(fff$names, collapse = ";")
              )
            }
          df.comb <- merge.data.frame(dfxx, df.names, by = "protein")
          df.comb$ontology <- ontology
          df.comb
        }
      colnames(mypathways)[2] <- "fold"
      mypathways$name.short <- paste(
        stringr::str_trunc(mypathways$protein, 15), stringr::str_trunc(mypathways$name, 40)
      )
      p <- strsplit(pp, "(", fixed = T)[[1]][1]
      ff <- data.frame(table(mypathways$ontology))
      ff <- ff[order(ff$Freq), ]
      plot.all <- ggplot(mypathways, aes(x = ontology, y = name.short)) +
        scale_color_manual(values = protools2::mycolors()$C23) +
        geom_point(aes(size = fold, color = ontology)) +
        theme_bw() +
        theme(axis.text.x = element_blank())
      list.ont <- list()
      ii <- 1
      for (ont in unique(mypathways$ontology)) {
        h <- mypathways[mypathways$ontology == ont, ]
        p <- ggplot(h, aes(y = reorder(name.short, fold), x = fold)) +
          geom_point(aes(color = -log(pvalues), size = -log(pvalues))) +
          labs(
            title = ont,
            color = "-log10(p)",
            size = "", y = ""
          ) +
          scale_color_gradient2(low = "seagreen", high = "purple4") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 1))
        list.ont[[ii]] <- list(data.table = h, dot.pot = p)
        names(list.ont)[[ii]] <- ont
        ii <- ii + 1
      }
      return(list(combined.plot = plot.all, list.of.plots.and.data = list.ont))
    }
    prots.in.increased.pathways <- .look.at.proteins.in.pathways(enrich.up)
    prots.in.decreased.pathways <- .look.at.proteins.in.pathways(enrich.do)
  }
  else {
    return(0)
  }
  return(list(
    pathway_enrichment_data = rbind.data.frame(enrich.do, enrich.up),
    delta_pathway_enrichment_data = rbind.data.frame(diff.enrich),
    plot_delta_enrichment = plot.diff,  # Used in report
    plot.pathways.enrichment.by.counts = plot.pathways.enrichment.by.counts,  # Used in report
    plot.pathways.enrichment.by.enrichment = plot.pathways.enrichment.by.enrichment,
    prots.in.increased.pathways = prots.in.increased.pathways,
    prots.in.decreased.pathways = prots.in.decreased.pathways)
  )
}


# Edited `protools2::plot.ont.vs.prot.relationships()` ----
# Modified to adjust y axis label padding to prevent overlapping of `cowplot` figure labels
plot.ont.vs.prot.relationships_edit <- function(results.ontology.analysis, n.pathways = 12, graph.title = "") {
  library(ggrepel)
  library(igraph)
  library(ggplot2)
  library(cowplot)
  .network.plot <- function() {
    pathways <- rresults$pathway[1:n.pathways]
    pth1 <- character()
    pth2 <- character()
    pval1 <- numeric()
    pval2 <- numeric()
    enrichment <- numeric()
    w <- numeric()
    g <- character()
    counts <- numeric()
    i <- 1
    df.network <- foreach(pathway = pathways, .combine = "rbind") %do% {
      count <- rresults[rresults$pathway == pathway, "counts"]
      pvalue <- rresults[rresults$pathway == pathway, "pvalues"]
      genes <- unlist(strsplit(rresults[rresults$pathway == pathway, "genes"], ";"))
      data.frame(
        pathway = rep(pathway, length(genes)),
        protein = genes
      )
    }
    g <- graph.data.frame(df.network, directed = T)
    V(g)$type <- bipartite.mapping(g)$type  # Occasionally produces error when `.network.plot()` called for "Decreased" (see below, line 458)
    # "Error in i_set_vertex_attr(x, attr(value, "name"), index = value, value = attr(value,  : Length of new attribute value must be 1 or 8, the number of target vertices, not 0"
    V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
    V(g)$shape <- ifelse(V(g)$type, "circle", "square")
    E(g)$color <- "lightgray"
    V(g)$label.color <- ifelse(V(g)$type, "black", "purple4")
    V(g)$label.cex <- ifelse(V(g)$type, 0.7, 1.2)
    V(g)$frame.color <- "gray"
    V(g)$size <- ifelse(V(g)$type, 3, 8)
    V(g)$size <- degree(g)
    fr.all <- layout.fruchterman.reingold(g)
    fr.all <- layout_with_graphopt(g)
    df <- data.frame(
      name = V(g)$name, color = V(g)$label.color,
      type = ifelse(V(g)$type, 1, 2), font.size = V(g)$label.cex,
      degree = degree(g, mode = "all")
    )
    fr.all.df <- as.data.frame(fr.all)
    fr.all.df$name <- V(g)$name
    gx <- get.data.frame(g)
    gx$from.x <- fr.all.df$V1[match(gx$from, fr.all.df$name)]
    gx$from.y <- fr.all.df$V2[match(gx$from, fr.all.df$name)]
    gx$to.x <- fr.all.df$V1[match(gx$to, fr.all.df$name)]
    gx$to.y <- fr.all.df$V2[match(gx$to, fr.all.df$name)]
    gx$name <- fr.all.df$name[match(gx$from, fr.all.df$name)]
    xxx <- merge.data.frame(fr.all.df, df, by = "name")
    pplot <- ggplot() +
      geom_segment(
        data = gx, aes(
          x = from.x, xend = to.x, y = from.y, yend = to.y
        ),
        colour = "grey4",
        linetype = 3, alpha = 0.5
      ) +
      geom_point(
        data = xxx, aes(
          x = V1, y = V2, shape = as.character(type),
          color = as.character(type), size = as.character(type)
        )
      ) +
      geom_text_repel(
        data = xxx, aes(
          x = V1, y = V2,
          label = name, colour = as.character(type), size = as.character(type)
        ),
        alpha = 1
      ) +
      scale_colour_manual(values = c(`2` = "purple4", `1` = "black")) +
      scale_size_manual(values = c(3, 5)) +
      scale_x_continuous(expand = c(0, 1)) +
      scale_y_continuous(expand = c(0, 1)) +
      labs(y = "") +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 30), # Modified axis.title.y
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), legend.position = "none"
      ) +
      ggtitle(my.title)
    pplot
    return(pplot)
  }
  x.inc <- subset(results.ontology.analysis, results.ontology.analysis$effect ==
                    "Increased")
  rresults <- x.inc[order(-x.inc$counts), ]
  my.title <- paste(graph.title, "Increased")
  plot.inc <- try (
    # plot.inc <-
    .network.plot(),
    silent = TRUE
  )
  x.dec <- subset(results.ontology.analysis, results.ontology.analysis$effect ==
                    "Decreased")
  rresults <- x.dec[order(-x.dec$counts), ]
  my.title <- paste(graph.title, "Decreased")
  plot.dec <- try (
    # plot.dec <-
    .network.plot(),  # Produces error (see comment above, line 384)
    silent = TRUE
  )
  # plots_list <- list(plot.inc, plot.dec)
  if (!inherits(plot.inc, "try-error") & !inherits(plot.dec, "try-error")) {
    return(list(
      plot.increased = plot.inc, plot.decreased = plot.dec,
      combined.plot = cowplot::plot_grid(plot.inc, plot.dec, nrow = 2)
    ))
  } else if (!inherits(plot.inc, "try-error") & inherits(plot.dec, "try-error")) {
    return(list(
      plot.increased = plot.inc,
      combined.plot = cowplot::plot_grid(plot.inc)
    ))
  } else if (inherits(plot.inc, "try-error") & !inherits(plot.dec, "try-error")) {
    return(list(
      plot.decreased = plot.dec,
      combined.plot = cowplot::plot_grid(plot.dec)
    ))
  } else {
    return(0)
  }
}


# Edited `protools2::kinase.substrate.enrichment()` ----
# To fix error:
# Error in `$<-.data.frame`(`*tmp*`, "prot.group", value = list(gene = "Akt1")) :
# replacement has 1 row, data has 9
kinase.substrate.enrichment_edit <- function(dfx, ks_db, is.ksea = TRUE) {
  dfx <- data.frame(dfx)
  # rownames(dfx) <- make.names(dfx[, 1], unique = T)  # No longer necessary
  if (is.ksea) {  # Original: `is.ksea == TRUE`
    column.with.priors <- 3
  }
  else {
    x1 <- grepl("(", as.character(dfx[, 1]), fixed = T)
    if (TRUE %in% x1) {
      dfx[, 1] <- unlist(lapply(
        dfx[, 1], function(x) strsplit(as.character(x), "(", fixed = TRUE)[[1]][1]
      ))
      column.with.priors <- "genes"
    }
  }
  # dfx[, 1] <- gsub("..", ");", dfx[, 1], fixed = T) # Unnecessary if `make.names()` not applied prior to initial `ksea()`
  # dfx[, 1] <- gsub(".", "(", dfx[, 1], fixed = T) # Unnecessary if `make.names()` not applied prior to initial `ksea()`
  nc <- ncol(dfx)
  df.ks <- protools2::protein_and_ks_sets[[ks_db]]
  nr <- nrow(df.ks)
  results.pvalues <- numeric(nr)
  results.zscores <- numeric(nr)
  results.q <- numeric(nr)
  results.m <- numeric(nr)
  results.distance <- numeric(nr)
  results.sites <- character(nr)
  kinases <- character(nr)
  r = 1
  for (r in 1:nr) {
    mym <- df.ks[r, 2]
    kinase <- df.ks[r, 1]
    if (!is.na(mym)) {  # Changed from `is.na(mym) == F` to `!is.na(mym)`
      if (mym >= 2) {  # Changed from original `> 2` to `>= 2`
        substrates <- as.character(df.ks[r, column.with.priors])
        ss <- c(unlist(strsplit(substrates, ";")))
        start.time <- Sys.time()
        # df.xx <- subset(dfx, dfx[, 1] %in% paste(ss, ";", sep = ""))  # Changed to `subset(dfx, dfx[, 1] %in% ss)` below
        df.xx <- subset(dfx, dfx[, 1] %in% ss)  # Changed from above
        sites <- paste(unlist(rownames(df.xx)), collapse = ";")
        sites <- gsub(";;", ";", sites) # Added
        if (nrow(df.xx) >= 2) {  # Changed from original `> 2` to `>= 2`
          df.xx$prot.group <- rep(kinase, nrow(df.xx))  # Edited to correct error, original: `<- kinase`
          sites.x <- data.frame(site = df.xx[, 1])
          sites.x$kinase <- rep(kinase, nrow(df.xx))  # Edited to correct error, original: `<- kinase`
          c = 2
          ds <- numeric(nc - 1)
          pvals <- numeric(nc - 1)
          zscores <- numeric(nc - 1)
          ms <- numeric(nc - 1)
          qs <- numeric(nc - 1)
          values.all <- na.omit(as.numeric(subset(dfx[, c], dfx[, c] != 0)))
          myvalues <- na.omit(as.numeric(subset(df.xx[, c], df.xx[, c] != 0)))
          pval <- 1
          tryCatch({
            myks <- ks.test(values.all, myvalues)
            pval <- myks$p.value
          },
          error = function(e) {}
          )
          m <- 0
          q <- 0
          m <- nrow(df.xx)
          mysd <- sd(values.all, na.rm = T)
          mymedian <- median(myvalues, na.rm = T)
          mymedian.all <- median(values.all, na.rm = T)
          sd.all <- sd(values.all)
          results.distance[r] <- mymedian - mymedian.all
          results.pvalues[r] <- pval
          results.zscores[r] <- ((mymedian - mymedian.all) * ((sqrt(m))) / mysd)
          results.m[r] <- m
          results.sites[r] <- as.character(sites)
          kinases[r] <- kinase
        }
      }
    }
  }
  xx <- data.frame(
    kinases, zscores = results.zscores, pvalues = results.pvalues,
    m = results.m, distance = results.distance, sites = results.sites,
    kinase_dataset = paste0(kinases, "_", ks_db)
  )
  xx$kinase_group <- unlist(kinases)
  xx$sites_kinase_dataset <- paste0(results.sites, "_", xx$kinase_dataset)  # To provide a column for unique row names
  xx <- subset(xx, xx$m > 1)
  xx$qvalue <- p.adjust(xx$pvalues, method = "fdr")
  return(xx)
}


# Edited `protools2::ksea()` ----
# To fix error in nested function:
# `protools2::kinase.substrate.enrichment()` (see above this function)
ksea_edit <- function(
    df.fold,
    ks_db = c("pdts", "psite"),
    graph.heading = "",
    pval.cut.off = 0.05
){
  library(ggrepel)
  yy <- protools2::expand.phosphopeptide.dataset(df.fold)
  cores = detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  xx <- foreach(db = ks_db, .combine = "rbind") %do% {
    # x <- protools2::kinase.substrate.enrichment(dfx = yy, ks_db = db)  # Edited, eventually build package to correct actual function
    x <- kinase.substrate.enrichment_edit(dfx = yy, ks_db = db)
    x
  }
  stopCluster(cl)
  # Labels
  if (ks_db == "ctams") {
    subtitle <- paste(
      "Changes in Compound Target Activity Markers (CTAMs)",
      "by Kinase Substrate Enrichment Analysis (KSEA)",
      sep = "\n"
    )
  } else {
    subtitle <- "Kinase Substrate Enrichment Analysis (KSEA)"
  }
  caption <-  paste(
    "Kinase activities quantified from changes in substrate phosphorylation.",
    "Labels indicate significant points.",
    sep = "\n"
  )
  df.k.up <- subset(xx, xx$pvalues < pval.cut.off & xx$zscores > 0)
  df.k.do <- subset(xx, xx$pvalues < pval.cut.off & xx$zscores < 0)
  plot.v.ksea.p <- ggplot(xx, aes(x = zscores, y = -log10(pvalues))) +
    geom_point() + geom_point(
      data = df.k.up,
      aes(x = zscores, y = -log10(pvalues), size = m),
      color = "red") +
    geom_label_repel(
      data = df.k.up,
      aes(x = zscores, y = -log10(pvalues), label = kinase_dataset),
      color = "red") +
    geom_point(
      data = df.k.do,
      aes(x = zscores, y = -log10(pvalues), size = m),
      color = "blue") +
    geom_label_repel(
      data = df.k.do,
      aes(x = zscores, y = -log10(pvalues), label = kinase_dataset),
      color = "blue") +
    labs(
      title = graph.heading,
      subtitle = subtitle,
      caption = caption) +
    geom_vline(xintercept = 0, linetype = 2) + theme_bw()
  df.k.up <- subset(xx, xx$qvalue < pval.cut.off & xx$zscores > 0)
  df.k.do <- subset(xx, xx$qvalue < pval.cut.off & xx$zscores < 0)
  plot.v.ksea.q <- ggplot(xx, aes(x = zscores, y = -log10(qvalue))) +
    geom_point() + geom_point(
      data = df.k.up,
      aes(x = zscores, y = -log10(qvalue), size = m),
      color = "red") +
    geom_label_repel(
      data = df.k.up,
      aes(x = zscores, y = -log10(qvalue), label = kinase_dataset),
      color = "red") +
    geom_point(
      data = df.k.do,
      aes(x = zscores, y = -log10(qvalue), size = m),
      color = "blue") +
    geom_label_repel(
      data = df.k.do,
      aes(x = zscores, y = -log10(qvalue), label = kinase_dataset),
      color = "blue") +
    labs(
      title = graph.heading,
      subtitle = subtitle,
      caption = caption) +
    geom_vline(xintercept = 0, linetype = 2) + theme_bw()
  return(
    list(
      ksea.data = xx,
      ksea.volcano.plot.pvals = plot.v.ksea.p,
      ksea.volcano.plot.qvals = plot.v.ksea.q
    )
  )
}

# Edited heatmap ----
#' Heatmap with p-value significance level labels.
#'
#' @param df.zcr A `data.frame` with fold changes or z-scores
#' @param df.pval A `data.frame` with p-values or q-values
#' @param mytitle A `string` for the heatmap title.
#' @param key.title A `string` for the key title, defaults to "log2Fold".
#'
#' @return A \code{\link{gplots::heatmap.2}} heatmap.
#' @export
myHeatMap_With_PValues <- function(
    df.zcr,
    df.pval,
    mytitle = "",
    key.title = "log2Fold"
) {

  # Colours
  palette.breaks <- seq(-1.5, 1.5, 0.1) # these values can be changed if needed
  color.palette <- colorRampPalette(c("blue", "white", "red"))(length(palette.breaks) - 1)

  # P-values
  df.pp <- data.frame(matrix(nrow = nrow(df.zcr), ncol = ncol(df.zcr)))
  for (cc in 1:ncol(df.pp)) {
    for (rr in 1:nrow(df.pp)) {
      vv <- df.pval[rr, cc]
      if (is.numeric(vv) & is.na(vv) == FALSE) {
        if (vv < 0.05) { df.pp[rr, cc] <- "*" }
        if (vv < 0.01) { df.pp[rr, cc] <- "**" }
        if (vv < 0.005) { df.pp[rr, cc] <- "***" }
        if (vv > 0.05) { df.pp[rr, cc] <- "" }
      } else {
        df.pp[rr, cc] <- ""
      }
    }
  }

  # Heatmap
  hm <- gplots::heatmap.2(
    as.matrix(df.zcr),
    cellnote = as.matrix(df.pp),
    notecex = 0.5,
    notecol = alpha("black", alpha = 0.5),
    dendrogram = "row",
    scale = "none",
    density.info = "none",
    symbreaks = FALSE,
    symkey = FALSE,
    Colv = NULL,
    # Rowv = NULL,
    trace = "none",
    margins = c(10, 20),
    col = color.palette,
    breaks = palette.breaks,
    cexCol = 1,
    cexRow = 0.5,
    colsep = 1:ncol(df.pp),
    rowsep = 1:nrow(df.pp),
    sepcolor = alpha("black", alpha = 0.5),
    sepwidth = c(0.01, 0.01),
    key = TRUE,
    key.title = NA,
    key.xlab = key.title,
    keysize = 0.8,
    # key.par=list(mar=c(1,0,1,1), cex=1.0, cex.lab=1.0, cex.axis=1.0),
    # lmat = rbind(c(0,3),c(2,1),c(0,3)),
    # lhei=c(0.1, 0.1),
    # lwid=c(1,4),
    # ColSideColors=df.samples,
    # ColSideColorsSize=10,
    main = mytitle
  )

  return(hm)
}


# ComplexHeatmap function ----
#' Heatmap with p-value significance level labels.
#'
#' @param zscores_df `data.frame`. Data containing fold changes or z-scores.
#' @param pvalues_df `data.frame`. Data containing p-values or q-values.
#' @param row_title `string`. The row title.
#' @param column_title `string`. The column title.
#' @param legend_title `string`. The key title.
#' @param fig_width `numeric`. The figure width (inches).
#' @param fig_height `numeric`. The figure height (inches). Default resizes
#'  dynamically according to number of rows.
#'
#' @return `Large Heatmap`. See \code{\link{ComplexHeatmap::Heatmap}}
#' @export
complex_heatmap <- function(
    zscores_df,
    pvalues_df,
    row_title = "",
    column_title = "",
    legend_title = "",
    fig_width=10,
    fig_height=NULL
) {

  # Calculate dynamic figure height
  if (is.null(fig_height)) {
    base_height <- 4  # base height when there are no rows
    height_per_row <- 0.1  # additional height per row
    fig_height <- base_height + (nrow(zscores_df) * height_per_row)
  }

  # Cell annotation
  df.pp <- data.frame(matrix(nrow = nrow(zscores_df), ncol = ncol(zscores_df)))
  for (cc in 1:ncol(df.pp)) {
    for (rr in 1:nrow(df.pp)) {
      vv <- pvalues_df[rr, cc]
      if (is.numeric(vv) & !is.na(vv)) {
        if (vv < 0.05) {
          df.pp[rr, cc] <- "*"
        }
        if (vv < 0.01) {
          df.pp[rr, cc] <- "**"
        }
        if (vv < 0.005) {
          df.pp[rr, cc] <- "***"
        }
        if (vv > 0.05) {
          df.pp[rr, cc] <- ""
        }
      } else {
        df.pp[rr, cc] <- ""
      }
    }
  }

  # Output format
  # pdf(
  #   file = "heatmap.pdf",
  #   width = fig_width, height = fig_height, compress = FALSE
  # )

  # Heatmap
  ComplexHeatmap::Heatmap(
    matrix = as.matrix(zscores_df),
    col = circlize::colorRamp2(
      breaks = c(
        min(zscores_df),
        0,
        max(zscores_df)
      ),
      # colors = c(
      #   viridis_pal(begin = 0, end = 1, option = "D")(3)[1], # hue_pal()(8)[8],
      #   viridis_pal(begin = 0, end = 1, option = "D")(3)[2], # "white",
      #   viridis_pal(begin = 0, end = 1, option = "D")(3)[3]  # hue_pal()(8)[4]
      # )
      colors = c(
        "blue",
        "white",
        "red"
      )
    ),
    name = legend_title,
    # Rows
    row_title = row_title,
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 9, fontface = "plain"),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 4),
    show_row_dend = TRUE,
    row_dend_side = "right",
    row_dend_width = unit(30, "mm"),
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    # Columns
    column_title = column_title,
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 9, fontface = "plain"),
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 7),
    show_column_dend = FALSE,
    # column_dend_side = "top",
    # column_dend_height = unit(30, "mm"),
    cluster_columns = FALSE,
    # clustering_distance_columns = "euclidean",
    # clustering_method_columns = "ward.D2",
    # Body
    # width =  fig_width - 0.5,   # Heatmap body
    # height =  fig_height - 1.5,  # Heatmap body
    heatmap_width = unit(fig_width, "inches"),  # Whole figure
    heatmap_height = unit(fig_height, "inches"),  # Whole figure
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      legend_direction = "horizontal", # "vertical",
      title_position = "topcenter",  # "leftcenter", # "leftcenter-rot",
      title = paste0("\n\n\n\n\n\n\n\n", legend_title),
      title_gp = gpar(fontsize = 5, fontface = "plain"),
      labels_gp = gpar(fontsize = 5)
    ),
    border = TRUE,
    border_gp = gpar(col = "black", alpha = 0.2),
    # Cell annotation
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        df.pp[i, j], x, y,
        gp = gpar(fontsize = 4)
      )
    }
  )
}


# UMAP function definition ----
#' Uniform Manifold Approximation and Projection (UMAP)
#'
#' Performs Uniform Manifold Approximation and Projection (UMAP) dimensionality
#' reduction on a matrix-like object using the \code{\link{umap::umap}} package.
#'
#' @param input \code{\link{data.frame}} or \code{\link{matrix}}-like object.
#' Contains numeric data.
#' @param label `'numeric'`, `'character'`, or `NULL`. Default `NULL` uses
#' \code{\link{rownames}}, else labels are taken from the specified column number
#' (`'numeric'`) or name (`'character'`) and separated from the numeric data.
#' @param size `'numeric'`. Label text size, default = `3`.
#' @param max.overlaps `numeric`. \code{\link{ggrepel::geom_label_repel}} `max.overlaps`
#' parameter to exclude text labels that overlap too many items; default = `10`.
#' @param title `'character'`. The title of the plot.
#' @param guide `bool`. Plots the legend if `TRUE`, default `FALSE` does not.
#'
#' @return \code{\link{list}}. A list of two elements:
#' \enumerate{
#'   \item \code{\link{umap::umap}} object. Containing at least a component with an
#'   embedding and a component with configuration settings.
#'   \item \code{\link{ggplot2::ggplot}} object. Plot of the UMAP result.
#' }
#' @export
#'
#' @section Notes:
#' The UMAP algorithm is stochastic, so repeated executions will yield different results.
#'
#' @section Warnings:
#' If output:
#' \preformatted{
#' `Warning: ggrepel: _ unlabelled data points (too many overlaps). Consider increasing max.overlaps
#' }
#' Resolve by any of the following:
#' \enumerate{
#'   \item Increase plot viewer size and/or increase chunk options `fig.width` and `fig.height`.
#'   \item Decrease argument of `size`.
#'   \item Increase argument of `max.overlaps`.
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{umap::umap}}
#' }
#'
#' @examples
umap_and_plot <- function(
    input,
    label = NULL,
    size = 3,
    max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
    title = "UMAP",
    guide = FALSE
) {
  # Packages
  require(umap)
  require(tidyverse)
  require(ggrepel)

  # Separate data
  matrix_data <- na.omit(input)
  if (is.null(label)) {
    matrix_labels <- rownames(matrix_data)
  } else if (is.numeric(label)) {
    matrix_labels <- matrix_data[, label]
    matrix_data <- matrix_data[, -label]
  } else if (is.character(label)) {
    matrix_labels <- matrix_data[, label]
    matrix_data <- matrix_data[, !(colnames(matrix_data) %in% label)]
  } else {
    stop(
      "Error: `label` must be one of `NULL`, of 'numeric' type, or of 'character' type.
      Hint: Specify a column in the matrix that contains label information. Default `NULL` uses `rownames`."
    )
  }

  # Projection
  message(paste("UMAP computation start:", format(Sys.time(), "%d-%m-%Y %H:%M:%S")))
  umap_projection <- umap::umap(
    matrix_data,
    n_neighbors = min(nrow(matrix_data) - 1, umap::umap.defaults$n_neighbors)
  )
  message(paste("UMAP computation end:", format(Sys.time(), "%d-%m-%Y %H:%M:%S")))

  umap_layout <- as.data.frame(umap_projection$layout)
  umap_layout <- umap_layout %>%
    dplyr::mutate(labels = matrix_labels) %>%
    dplyr::rename(UMAP1 = V1, UMAP2 = V2)

  # Plot
  plot_umap <- ggplot(
    umap_layout,
    aes(
      x = UMAP1,
      y = UMAP2,
      colour = labels
    )
  ) +
    geom_point() +
    geom_label_repel(
      aes(label = labels),
      size = size,
      max.overlaps = max.overlaps
    ) +
    theme_bw() +
    labs(
      title = title
    )
  # Omit colour legend
  if (!guide) {
    plot_umap <- plot_umap + guides(colour = "none")
  }
  return(
    list(umap_projection = umap_projection, plot_umap = plot_umap)
  )
}
