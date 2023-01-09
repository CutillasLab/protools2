


plot.kinase.relationships <- function(ksea.data, graph.title, pval.cutoff){

  library(ggrepel)
  library(igraph)
  kinases <- ksea.data[ksea.data$pvalues<pval.cutoff,]$kinase_dataset
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
  for (k1 in kinases){

    count <- ksea.data[ksea.data$kinase_dataset==k1,"m"]
    pvalue <- ksea.data[ksea.data$kinase_dataset==k1,"pvalues"]
    zscore <- ksea.data[ksea.data$kinase_dataset==k1,"zscores"]
    #if (count>1 & pvalue<0.05){
    genes.1 <- unlist(strsplit( ksea.data[ksea.data$kinase_dataset==k1,"sites"],";"))
    for (k2 in kinases){
      genes.2 <- unlist(strsplit( ksea.data[ksea.data$kinase_dataset==k2,"sites"],";"))

      #if (pathway.1!=pathway.2){
      common <- intersect(genes.1,genes.2)
      if (length(common)>0){
        pth1[i] <- k1
        pth2[i] <- k2
        g[i] <- paste0( common,collapse = ";")
        w[i] <- length(common)
        pval1[i] <-  ksea.data[ksea.data$kinase_dataset==k2,"pvalues"]
        pval2[i] <-  ksea.data[ksea.data$kinase_dataset==k1,"pvalues"]
        zscores[i] <- ksea.data[ksea.data$kinase_dataset==k1,"zscores"]
        counts[i] <- ksea.data[ksea.data$kinase_dataset==k1,"m"]
        i <- i+1
      }
      #}
      #}
    }

  }

  df.network <- data.frame(node.1 =pth1,
                           node.2=pth2,
                           weight=w,
                           #genes=g,
                           pvalue.1 =pval1,
                           pvalue.2=pval2,
                           zscores=zscores,
                           counts=counts)

  df.network <- df.network[order(df.network$pvalue.1),]
  df.s <- subset(df.network,df.network$weight>4)

  xx <- reshape2::dcast(df.s, node.1~node.2, value.var = 'weight')
  xx[is.na(xx)] <- 0
  rownames(xx) <- xx$node.1

  caught.inc <- graph.incidence(xx[,2:ncol(xx)], weighted = TRUE)  #make data into a bipartite graph object
  obs.parties.all <- bipartite.projection(caught.inc)[[1]]
  obs.spp.all <- bipartite.projection(caught.inc)[[2]]

  fr.all <- layout_with_graphopt(obs.spp.all)

  fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
  fr.all.df$species <- colnames(xx[,2:ncol(xx)])  ## add in the species codes

  g <- get.data.frame(obs.spp.all)
  g$from.x <- fr.all.df$V1[match(g$from, fr.all.df$species)]  #  match the from locations from the node data.frame we previously connected
  g$from.y <- fr.all.df$V2[match(g$from, fr.all.df$species)]
  g$to.x <- fr.all.df$V1[match(g$to, fr.all.df$species)]  #  match the to locations from the node data.frame we previously connected
  g$to.y <- fr.all.df$V2[match(g$to, fr.all.df$species)]

  g$pvalue.1 <- df.network$pvalue.1[match(g$from,df.network$node.1)]
  g$pvalue.2 <- df.network$pvalue.2[match(g$from,df.network$node.2)]
  g$weight <- df.network$weight[match(g$from,df.network$node.1)]


  fr.all.df$pvalue.1 <- df.network$pvalue.1[match(fr.all.df$species,df.network$node.1)]
  fr.all.df$pvalue.2 <- df.network$pvalue.2[match(fr.all.df$species,df.network$node.2)]



  fr.all.df$E <- df.network$zscores[match(fr.all.df$species,df.network$node.1)]
  fr.all.df$counts <- df.network$counts[match(fr.all.df$species,df.network$node.2)]

  pplot <- ggplot() +
    geom_segment(data=g,
                 aes(x=from.x,xend = to.x, y=from.y,yend = to.y, size=weight/10),colour="black",linetype=2) + # add line type
    #geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes
    geom_point(data=fr.all.df,aes(x=V1,y=V2,fill=E,size=counts),shape=21) +
    geom_label_repel(data=fr.all.df,aes(x=V1,y=V2,label=species)) + # add the node labels
    #scale_colour_manual(values=c("1"="red","2"="lightblue"))+  # add colour scaling for group membership
    scale_fill_gradient2(low="blue",high = "red",mid = "white",midpoint = 0)+
    #scale_linetype_manual(values=c("0"="dashed","1"="solid"))+ # add linteyp scaling for within and between groups
    scale_x_continuous(expand=c(0,1))+  # expand the x limits
    scale_y_continuous(expand=c(0,1))+ # expand the y limits
    theme_bw()+  # use the ggplot black and white theme
    theme(
      axis.text.x = element_blank(),  # remove x-axis text
      axis.text.y = element_blank(), # remove y-axis text
      axis.ticks = element_blank(),  # remove axis ticks
      axis.title.x = element_blank(), # remove x-axis labels
      axis.title.y = element_blank(), # remove y-axis labels
      panel.background = element_blank(),
      panel.border =element_blank(),
      panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),  #remove minor-grid labels
      plot.background = element_blank())+
    ggtitle(graph.title)+
    labs(caption = "Edge weights are proportional to common downstream targets between kinases. \nNode colors are proportional to z-score of enrichment.")
  pplot
  return(pplot)

}
