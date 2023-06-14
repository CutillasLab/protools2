# Functions to facillitate Proteomics and Phosphoproteomics script execution


# Fix incorrect combiPeptData gene names ----
#' Fixes incorrect data of pre-processed mass spectrometry input:
#' combiPeptData columns 25 and 29.
#'
#' @param df.combi \code{\link{data.frame}}. combiPeptData.
#' @param organism String. One of "human" or "mouse". Selects database of
#' Uniprot names.
#' @param pescal.xlsx String. The path ("directory/name") of the Pescal Excel
#' output.
#' @param fixed.xlsx String. Name for the fixed output Excel file.
#'
#' @return \code{\link{data.frame}}. The corrected combiPeptData.
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
    df.uniprot <- protools2::uniprot.names
  } else if (organism == "mouse") {
    df.uniprot <- protools2::uniprot.names.mouse
  } else {
    stop("Error: `organism` must be one of 'human' or 'mouse'")
  }

  # Fix genes column
  df1 <- df.combi %>%
    mutate(uniprot_split = str_split(uniprot, "; ")) %>%
    unnest(uniprot_split) %>%
    mutate(uniprot_split = str_replace_all(uniprot_split, ";", ""))

  df2 <- df.uniprot %>%
    filter(Entry %in% df1$uniprot_split) %>%
    select(Entry, Gene.Names..primary.)

  df3 <- df1 %>%
    left_join(df2, by = c("uniprot_split" = "Entry")) %>%
    group_by(db_id) %>%
    summarise(genes_new = paste(Gene.Names..primary., collapse = "; ")) %>%
    mutate(genes_new = paste0(genes_new, ";"))

  df4 <- df.combi %>%
    left_join(df3, by = "db_id")

  # Fix genes_mod column
  df5 <- df4 %>%
    rowwise() %>%
    mutate(
      genes_mod_vals = list(str_split_1(genes_mod, ";")),
      genes_mod_keep = list(unlist(str_extract_all(genes_mod_vals, "\\(.*?\\)"))),
      genes_new_split = list(str_replace_all(str_split_1(genes_new, "; "), ";", "")),
      genes_mod_reinsert = list(paste0(genes_new_split, genes_mod_keep)),
      genes_mod_recomb = paste0(genes_mod_reinsert, collapse = ";"),
      genes_mod_recomb = paste0(genes_mod_recomb, ";")
    ) %>%
    mutate(
      genes = genes_new,
      genes_mod = genes_mod_recomb,
      genes_new = NULL,
      genes_mod_vals = NULL,
      genes_mod_keep = NULL,
      genes_new_split = NULL,
      genes_mod_reinsert = NULL,
      genes_mod_recomb = NULL
    ) %>%
    rename(
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
#' X <- data.table(x=c("c","b"), v=8:7, foo=c(4,2))
#'
#' multi_DT <- mergeDTs(
#'   list(DT, X),
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
