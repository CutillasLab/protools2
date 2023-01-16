


compare.conditions.by.ttest <- function(df.design,
                                        df.comparisons,
                                        df.normalized.areas,
                                        is.paired=FALSE)
  {


  df.comparisons$comparison <- paste0(df.comparisons$numerator,"_vs_",df.comparisons$denominator)
  comparisons <- df.comparisons$comparison
  comp <- comparisons[1]
  list.of.comparisons <- list()
  i <- 1
  for ( comp in comparisons){

    control.condi <- df.comparisons[df.comparisons$comparison==comp,"denominator"]$denominator
    test.condi <- df.comparisons[df.comparisons$comparison==comp,"numerator"]$numerator

    control.samples <- df.design[df.design$condition==control.condi,"heading"]$heading
    test.samples <- df.design[df.design$condition==test.condi,"heading"]$heading

    x <- compare.two.conditions.by.ttest(df.to.compare = df.normalized.areas,
                                         control.columns = control.samples,
                                         test.columns = test.samples,
                                         is.paired = is.paired)

    list.of.comparisons[[i]] <- x
    names(list.of.comparisons)[i] <- comp
    i <- i+1
  }
  return(list.of.comparisons)

}


compare.two.conditions.by.ttest <- function(df.to.compare,
         log.data=FALSE,
         control.columns="",
         test.columns="",
         is.paired=FALSE)
  {


  # Compare by ttest

  library(foreach)
  library(doParallel)

  nc <- ncol(df.to.compare)
  nr <- nrow(df.to.compare)

  if (log.data){
    df.s <- data.frame(prot=df.to.compare[,1], scale(log2(df.to.compare[,2:nc])))
  }else{
    df.s <- df.to.compare
  }

  if (length(control.columns)<2){
    print("Give me the names of control.columns and test.columns.")
    return()
  }
  if (length(test.columns)<2){
    print("Give me the names of control.columns and test.columns.")
    return()
  }

  registerDoParallel(cores = parallel::detectCores()-1)
  # print(paste("running on", parallel::detectCores()-1, "cores"))
  t1 <- Sys.time()

  modeled.values <- foreach(r = 1:nr, .combine = "rbind")%dopar%
    {

      cont <- unlist(df.s[r,control.columns])
      test <- unlist(df.s[r,test.columns])



      if (length(na.omit(cont))>2 & length(na.omit(test))>2){
        pval <- 1
        tryCatch(
        {pval <- t.test( cont,test, paired=is.paired )$p.value},error=function(e){}
        )
        fold <- mean(test,na.rm = T)-mean(cont,na.rm = T)
        result <- c(df.s[r,1],fold,pval)
      }else{
        result <- c(df.s[r,1],0,1)
      }
      result
    }
  t2 <- Sys.time()
  # print(t2-t1)
  stopImplicitCluster()

  df <- data.frame(modeled.values)
  colnames(df) <- c("protein","difference.test.vs.control","pvalues")
  df$difference.test.vs.control <- as.numeric(df$fold)
  df$pvalues <- as.numeric(df$pvalues)
  df$FDR <- p.adjust(df$pvalue,method = "fdr")
  rownames(df) <- make.names(df$protein,unique = T)
  return(df)
}
