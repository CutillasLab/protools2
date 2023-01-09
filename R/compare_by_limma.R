
compare.conditions.by.limma <- function(df.design, df.comparisons, df.normalized.areas){


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

    x <- compare.by.limma(df.normalized.areas,control.samples = control.samples, test.samples = test.samples)

    list.of.comparisons[[i]] <- x
    names(list.of.comparisons)[i] <- comp
    i <- i+1
  }
  return(list.of.comparisons)

}





compare.by.limma <- function(df.to.compare, control.samples, test.samples){


  # Compare by limma
  #
  # first column is protein or ppsite names
  # first set of samples are control
  # second set of samples are the test
  library(limma)
  # df.to.compare <- df.ppindex[,1:9]
  nc <- ncol(df.to.compare)
  #replicates <- (nc-1)/2
  control.samples <- intersect(control.samples,colnames(df.to.compare))
  test.samples <- intersect(test.samples,colnames(df.to.compare))
  df.s <- df.to.compare[,c(control.samples,test.samples)]
  df.s1 <- data.frame(outcome=matrix(nrow=length(control.samples)))
  df.s2 <- data.frame(outcome=matrix(nrow=length(test.samples)))
  df.s1$outcome <- "control"
  df.s2$outcome <- "test"
  df.ss <- rbind(df.s1,df.s2)
  des <- factor(ifelse(df.ss$outcome=="control" ,"1",
                       "2"))
  facna <- addNA(des)
  design <- model.matrix(~ 0+factor(c(facna)))
  colnames(design) <- c("control","test")
  contrast.matrix <- makeContrasts(test-control,
                                   levels=design)
  #df.s ==== proteomics data
  fit <- lmFit(df.s,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  pvals <- data.frame(fit2$p.value)
  fvals <- data.frame(fit2$coefficients)
  df.xx <- data.frame(protein=df.to.compare[,1],
                      differnces=fvals,
                      pvalues=pvals)
  colnames(df.xx) <- c("protein","difference.test.vs.control","pvalues")

  df.xx$FDR <- p.adjust(df.xx$pvalues, method = "fdr")

  return(df.xx)

}
