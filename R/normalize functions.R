# Normalise output_areas
# Author & copyright: Pedro R. Cutillas


# Impute missing values ----
# NOTE: Function uses vectorised operations for speed, loops only on columns not rows.
impute_na <- function(df.design, df.areas) {
  # Get unique conditions
  conditions <- unique(df.design$condition)

  # Loop over each unique condition
  for(condition in conditions) {
    # Get the columns in df.areas that belong to the current condition
    cols <- df.design$heading[df.design$condition == condition]

    # Calculate the mean of each row for the current condition, ignoring NA values
    row_means <- apply(df.areas[, cols], 1, function(x) if(all(is.na(x))) NA else mean(x, na.rm = TRUE))

    # Replace NA values in df.areas for the current condition with the row means
    for(heading in cols) {
      is_na <- is.na(df.areas[, heading])
      df.areas[is_na, heading] <- row_means[is_na]

      # If all values in the row for the current condition are NA, replace NA with the minimum value of the column - 1
      all_na <- is.na(row_means)
      df.areas[all_na, heading] <- min(df.areas[, heading], na.rm = TRUE) - 1
    }
  }

  return(df.areas)
}


# Normalise phospho ----
normalize_areas_return_ppindex <- function(
    pescal_output_file,
    delta_score_cut_off = 0  # if (fragpipe) {0} else {5}
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
  suppressMessages(
    df.design <- readxl::read_excel(pescal_output_file, "design")
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

  # Alternative NA imputation
  # Centred + scaled
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

  # New improved NA imputation
  df.norm.log2.centered.scaled.na.imputed.new <- impute_na(
    df.design = df.design,
    df.areas = df.norm.log2.centered.scaled
  )

  # Centred
  df.norm.log2.centered.na.imputed <- df.norm.log2.centered

  df.norm.log2.centered.na.imputed[cols] <- lapply(
    df.norm.log2.centered.na.imputed[cols], function(x){
      replace(x, is.na(x), min(x, na.rm = TRUE) -1) # Correct NA imputation
    }
  )

  # New improved NA imputation
  df.norm.log2.centered.na.imputed.new <- impute_na(
    df.design = df.design,
    df.areas = df.norm.log2.centered
  )

  # Previous NA imputation
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
    df.norm.log2.centered.na.imputed = df.norm.log2.centered.na.imputed,
    df.norm.log2.centered.na.imputed.new = df.norm.log2.centered.na.imputed.new,
    df.norm.log2.centered.scaled.na.imputed = df.norm.log2.centered.scaled.na.imputed,
    df.norm.log2.centered.scaled.na.imputed1 = df.norm.log2.centered.scaled.na.imputed1,
    df.norm.log2.centered.scaled.na.imputed2 = df.norm.log2.centered.scaled.na.imputed2,
    df.norm.log2.centered.scaled.na.imputed.new = df.norm.log2.centered.scaled.na.imputed.new
  ))
}


# Normalise total ----
normalize_areas_return_protein_groups <- function(
    pescal_output_file,
    mascot.score.cut.off = 50,  # if (fragpipe) {0} else {40}
    n.peptide.cut.off = 1
) {
  suppressMessages(
    df.areas <- data.frame(readxl::read_excel(pescal_output_file, "output_areas"))
  )
  colnames(df.areas) <- gsub("-", ".", colnames(df.areas), fixed = T)
  suppressMessages(
    df.combi <- data.frame(readxl::read_excel(pescal_output_file, "combiPeptData"))  # Original `pescal.output.file` object is accessed from global scope
  )
  suppressMessages(
    df.design <- readxl::read_excel(pescal_output_file, "design")
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

  # Centred
  # New improved NA imputation
  df.norm.log2.centered.na.imputed.new <- impute_na(
    df.design = df.design,
    df.areas = df.norm.log2.centered
  )

  # Centred + scaled
  df.norm.log2.centered.scaled.na.imputed <- df.norm.log2.centered.scaled

  # df.norm.log2.centered.scaled.na.imputed[  # Previous na imputation
  #   is.na(df.norm.log2.centered.scaled.na.imputed)
  # ] <- min(df.norm.log2.centered.scaled.na.imputed[, cols], na.rm = T) / 5

  # Correct NA imputation
  df.norm.log2.centered.scaled.na.imputed[cols] <- lapply(
    df.norm.log2.centered.scaled.na.imputed[cols], function(x){
      replace(x, is.na(x), min(x, na.rm = TRUE) -1) # Correct NA imputation
    }
  )

  # New improved NA imputation
  df.norm.log2.centered.scaled.na.imputed.new <- impute_na(
    df.design = df.design,
    df.areas = df.norm.log2.centered.scaled
  )

  rownames(df.norm) <- df.norm.log2.centered$protein.group
  rownames(df.norm.log2.centered) <- df.norm.log2.centered$protein.group
  rownames(df.norm.log2.centered.scaled) <- df.norm.log2.centered.scaled$protein.group
  rownames(df.norm.log2.centered.na.imputed.new) <- df.norm.log2.centered.na.imputed.new$protein.group
  rownames(df.norm.log2.centered.scaled.na.imputed) <- df.norm.log2.centered.scaled.na.imputed$protein.group
  rownames(df.norm.log2.centered.scaled.na.imputed.new) <- df.norm.log2.centered.scaled.na.imputed.new$protein.group

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
  df.norm.log2.centered.na.imputed.new <- merge.data.frame(
    df.norm.log2.centered.na.imputed.new,
    x[, cc],
    by = "protein.group"
  )
  df.norm.log2.centered.scaled.na.imputed.new <- merge.data.frame(
    df.norm.log2.centered.scaled.na.imputed.new,
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
      df.norm.log2.centered.scaled.na.imputed$protein.group %in% selected.prot.groups, ],
    df.norm.log2.centered.na.imputed.new = df.norm.log2.centered.na.imputed.new[
      df.norm.log2.centered.na.imputed.new$protein.group %in% selected.prot.groups, ],
    df.norm.log2.centered.scaled.na.imputed.new = df.norm.log2.centered.scaled.na.imputed.new[
      df.norm.log2.centered.scaled.na.imputed.new$protein.group %in% selected.prot.groups, ]
  ))
}


# Means t-test ----
normalize_df_to_row_mean_and_get_one_sample_ttest <- function(df){

  df.s <- data.frame(t(apply(df,1,function(x) x-median(x,na.rm=T))))

  .one.sample.ttest <- function(j){

    pvalues <- numeric()
    i <- 1
    # x <- NA
    for (x in j){
      if (is.na(x)==F){
        pvalues[i] <-   t.test(j, mu = x, alternative = "two.sided")$p.value
      }else{
        pvalues[i] <- NA
      }
      i <- i+1
    }
    return(pvalues)
  }

  yy <- data.frame(t(apply(df.s,1, function(v) one.sample.ttest(v) )))

  rownames(yy) <- rownames(df)

  return(list(
    df_normalize_to_meanRow=df.s,
    df_pvalues=yy
  ))
}

