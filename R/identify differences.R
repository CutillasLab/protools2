

identify_differences_in_comparison_plus_volcano <- function(df, fold.cutoff=0.5,qval.cutoff=0.05, graph.header=""){

  if (ncol(df)==4){
    colnames(df) <- c("sites","fold","pvalue","qvalue")
  }else{
    colnames(df) <- c("sites","fold","pvalue","qvalue", "accession")
  }

  df.up <- subset(df, df$qvalue<qval.cutoff & df$fold>fold.cutoff)
  df.do <- subset(df, df$qvalue<qval.cutoff & df$fold<(-fold.cutoff))
  df.up <- df.up[order(-df.up$fold),]
  df.do <- df.do[order(df.do$fold),]
  suppressMessages(
  plot.pval.dist <- ggplot(df,aes(x=pvalue))+
    geom_histogram(fill="white",color="black")+
    labs(title = "Distrubtion of pvalues")
  )
  suppressMessages(
  plot.qval.dist <- ggplot(df,aes(x=qvalue))+
    geom_histogram(fill="white",color="black") +
    labs(title = "Distrubtion of qvalues (BH adjusted pvalues)")
  )
  #### Volcano plot

  plot.volcano.1 <- ggplot(df,aes(x=fold,y=-log10(qvalue)))+
    geom_point(size=1)+
    geom_point(data=df.up,aes(x=fold,y=-log10(qvalue)),color="red")+
    geom_point(data=df.do,aes(x=fold,y=-log10(qvalue)),color="blue")+
    labs(title = graph.header, subtitle = paste0("Increased = ", nrow(df.up), "; Decreased = ", nrow(df.do)))

  return(list( df.increased=df.up,
               df.decreased=df.do,
               volcanoplot=plot.volcano.1,
               pvalue.distributions=cowplot::plot_grid(plot.pval.dist,plot.qval.dist)
               )
         )


}


identify_differences_by_pvalue_in_comparison_plus_volcano <- function(df, fold.cutoff=0.5,pval.cutoff=0.05, graph.header=""){

  if (ncol(df)==4){
  colnames(df) <- c("sites","fold","pvalue","qvalue")
  }else{
    colnames(df) <- c("sites","fold","pvalue","qvalue", "accession")
  }

  df.up <- subset(df, df$pvalue<pval.cutoff & df$fold>fold.cutoff)
  df.do <- subset(df, df$pvalue<pval.cutoff & df$fold<(-fold.cutoff))
  df.up <- df.up[order(-df.up$fold),]
  df.do <- df.do[order(df.do$fold),]
  suppressMessages(
    plot.pval.dist <- ggplot(df,aes(x=pvalue))+
      geom_histogram(fill="white",color="black")+
      labs(title = "Distrubtion of pvalues")
  )
  suppressMessages(
    plot.pval.dist <- ggplot(df,aes(x=pvalue))+
      geom_histogram(fill="white",color="black") +
      labs(title = "Distrubtion of pvalues (BH adjusted pvalues)")
  )
  #### Volcano plot

  plot.volcano.1 <- ggplot(df,aes(x=fold,y=-log10(pvalue)))+
    geom_point(size=1)+
    geom_point(data=df.up,aes(x=fold,y=-log10(pvalue)),color="red")+
    geom_point(data=df.do,aes(x=fold,y=-log10(pvalue)),color="blue")+
    labs(title = graph.header, subtitle = paste0("Increased = ", nrow(df.up), "; Decreased = ", nrow(df.do)))

  return(list( df.increased=df.up,
               df.decreased=df.do,
               volcanoplot=plot.volcano.1,
               pvalue.distributions=cowplot::plot_grid(plot.pval.dist,plot.pval.dist)
  )
  )


}
