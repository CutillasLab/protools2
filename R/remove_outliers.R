

remove_outlier_samples_from_dataset <- function(df){
  cols <- colnames(dplyr::select_if(df, is.numeric))
  cms <- apply(df[,cols],2,median)
  OutVals = boxplot(cms)$out
  samples.to.be.removed <- names(OutVals)
  df <- df[ , -which(names(df) %in% samples.to.be.removed)]
  return(df)
}

