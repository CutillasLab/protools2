





make.plots.of.significant.peptides <- function(list.of.comparisons,
                                               proteomic.data,
                                               df.comparisons,
                                               df.design,
                                               fold_cutoff,
                                               qval_cutoff){


  list.of.comparisons <- comparisons


  df.all <- do.call("rbind",list.of.comparisons)

    suppressMessages(
      ddd <- identify_differences_in_comparison_plus_volcano(df.all,
                                                             fold.cutoff = fold_cutoff,
                                                             qval.cutoff = qval_cutoff,
                                                             graph.header=names(dd))
    )
    length(unique(ddd$df.increased[,1]))
    length(unique( ddd$df.decreased[,1]))
    length(unique( c(ddd$df.decreased[,1], ddd$df.increased[,1])))

}

make.plots.of.peptides <- function(list.of.peptides,
                                               proteomic.data,
                                               df.comparisons,
                                               df.design){


  comps.to.make <- list()
  i <- 1
  for(r in 1:nrow(df.comparisons)){
    n <- df.comparisons$numerator[r]
    d <- df.comparisons$denominator[r]
    comps.to.make[[i]] <- c(n,d)
    i <- i+1
  }

  list.of.plots <- list()

  pep <- "AKT1S1(S247)"
  i <- 1
  for (pep in list.of.peptides){
    dfx <- data.frame(t(df.norm[pep,]))
    dfx <- merge.data.frame(dfx,df.design, by.x=0,by.y = "heading")
    colnames(dfx)[2] <- "ppindex"
    pepplot <- ggplot(dfx,aes(x=condition,y=as.numeric(ppindex)))+
      geom_boxplot()+
      theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      stat_compare_means(comparisons = comps.to.make)+
      labs(title = pep, y="Abundance (a.u.)")+
      ylim(c(0,NA))
    list.of.plots[[i]] <- pepplot
    names(list.of.plots)[i] <- pep
    i <- i+1
  }
  return(list.of.plots)
}

