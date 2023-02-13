

get.ks.relationships.graph <- function(){

  # No not use , under development

ksea.data <- ke$ksea.data
ksea.data$alpha <- ksea.data$zscores*-log10(ksea.data$pvalues)
ksea.data.up <- ksea.data[order(-ksea.data$alpha),]
ksea.data.do <- ksea.data[order(ksea.data$alpha),]
kinases <- ksea.data[ksea.data$pvalues<0.05,]$kinase_dataset
k <- kinases[1]
df.network <- foreach(k = kinases,.combine = "rbind")%do%{

  #if (count>1 & pvalue<0.05){

  xx <- ksea.data[ksea.data$kinase_dataset==k,]
  genes <- unlist(strsplit(xx$sites,";"))

  data.frame(kinase=rep(k,length(genes)),
             protein=genes,
             zscore=xx$zscores,
             pvalue=xx$pvalues
             )


}

g <- graph.data.frame(df.network, directed = T)

V(g)$name
V(g)$type <- bipartite.mapping(g)$type
#V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
#V(g)$shape <- ifelse(V(g)$type, "circle", "square")
#E(g)$color <- "lightgray"
V(g)$label.color <- ifelse(V(g)$type, "black", "purple4")
## V(g)$label.font <-  2
#V(g)$label.cex <- ifelse(V(g)$type, 0.7, 1.2)
## V(g)$label.dist <-0
#V(g)$frame.color <-  "gray"
#V(g)$size <- ifelse(V(g)$type,3, 8)
#V(g)$size <- degree(g)

fr.all <- layout_with_graphopt(g)


df <- data.frame(name= V(g)$name,
                 color=V(g)$label.color,
                 type=ifelse(V(g)$type,1,2)
                 #,
                 #font.size=V(g)$label.cex,
                 #degree=degree(g, mode="all")
)

fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
fr.all.df$name <- V(g)$name# c(df.network$pathway,df.network$protein)# colnames(xx[,2:ncol(xx)])  ## add in the species codes

gx <- get.data.frame(g)
gx$from.x <- fr.all.df$V1[match(gx$from, fr.all.df$name)]  #  match the from locations from the node data.frame we previously connected
gx$from.y <- fr.all.df$V2[match(gx$from, fr.all.df$name)]
gx$to.x <- fr.all.df$V1[match(gx$to, fr.all.df$name)]  #  match the to locations from the node data.frame we previously connected
gx$to.y <- fr.all.df$V2[match(gx$to, fr.all.df$name)]
gx$name <- fr.all.df$name[match(gx$from, fr.all.df$name)]
xxx <- merge.data.frame(fr.all.df ,df, by="name")
xxx <- merge.data.frame(xxx, gx,by.x = "name",by.y = "from")

pplot <- ggplot() +
  geom_segment(data=gx,
               aes(x=from.x,xend = to.x, y=from.y,yend = to.y),
               colour="black",linetype=2) + # add line type
  #geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes

  geom_point(data=xxx,aes(x=V1,y=V2,shape=as.character(type),
                          color=zscore, size=zscore)) +


  geom_label_repel(data=xxx,aes(x=V1,y=V2,label=ifelse(type==2, name.y,"")),
                               colour="purple4",size=3) + # add the node labels
  #scale_colour_manual(values=c("2"="purple4","1"="black"))+  # add colour scaling for group membership
  #scale_size_manual(values=c(0.5,5))+
  #scale_color_gradient2(low="blue",high = "purple",mid = "white")+
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
    plot.background = element_blank(),
    legend.position = "none")+
  ggtitle(my.title)

plot(pplot)

}

plot.ont.vs.prot.relationships <- function(results.ontology.analysis,
                                           n.pathways=12,
                                           graph.title=""){

  library(ggrepel)
  library(igraph)
  library(ggplot2)
  library(cowplot)
  .network.plot <- function(){

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
    df.network <- foreach(pathway = pathways,.combine = "rbind")%do%{

      count <- rresults[rresults$pathway==pathway,"counts"]
      pvalue <- rresults[rresults$pathway==pathway,"pvalues"]
      #if (count>1 & pvalue<0.05){
        genes <- unlist(strsplit( rresults[rresults$pathway==pathway,"genes"],";"))
      data.frame(pathway=rep(pathway,length(genes)),protein=genes)
    }
    g <- graph.data.frame(df.network, directed = T)
    V(g)$type <- bipartite.mapping(g)$type
    V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
    V(g)$shape <- ifelse(V(g)$type, "circle", "square")
    E(g)$color <- "lightgray"
    V(g)$label.color <- ifelse(V(g)$type, "black", "purple4")
    ## V(g)$label.font <-  2
    V(g)$label.cex <- ifelse(V(g)$type, 0.7, 1.2)
    ## V(g)$label.dist <-0
    V(g)$frame.color <-  "gray"
    V(g)$size <- ifelse(V(g)$type,3, 8)
    V(g)$size <- degree(g)
    fr.all <- layout.fruchterman.reingold(g)
    fr.all <- layout_with_graphopt(g)
    df <- data.frame(name= V(g)$name,
                    color=V(g)$label.color,
                    type=ifelse(V(g)$type,1,2),
                    font.size=V(g)$label.cex,
                    degree=degree(g, mode="all")
    )

    fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
    fr.all.df$name <- V(g)$name# c(df.network$pathway,df.network$protein)# colnames(xx[,2:ncol(xx)])  ## add in the species codes

    gx <- get.data.frame(g)
    gx$from.x <- fr.all.df$V1[match(gx$from, fr.all.df$name)]  #  match the from locations from the node data.frame we previously connected
    gx$from.y <- fr.all.df$V2[match(gx$from, fr.all.df$name)]
    gx$to.x <- fr.all.df$V1[match(gx$to, fr.all.df$name)]  #  match the to locations from the node data.frame we previously connected
    gx$to.y <- fr.all.df$V2[match(gx$to, fr.all.df$name)]
    gx$name <- fr.all.df$name[match(gx$from, fr.all.df$name)]
    xxx <- merge.data.frame(fr.all.df ,df, by="name")
  pplot <- ggplot() +
      geom_segment(data=gx,
                   aes(x=from.x,xend = to.x, y=from.y,yend = to.y),
                   colour="grey4",linetype=3, alpha=0.5) + # add line type
      #geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes
      geom_point(data=xxx,aes(x=V1,y=V2,shape=as.character(type),
                              color=as.character(type),size=as.character(type))) +
      geom_text_repel(data=xxx,aes(x=V1,y=V2,label=name,
                                   colour=as.character(type),size=as.character(type)),alpha=1) + # add the node labels
      scale_colour_manual(values=c("2"="purple4","1"="black"))+  # add colour scaling for group membership
      scale_size_manual(values=c(3,5))+
  #scale_color_gradient2(low="blue",high = "purple",mid = "white")+
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
        plot.background = element_blank(),
        legend.position = "none")+
    ggtitle(my.title)
pplot
    return(pplot)
  }

  x.inc <- subset(results.ontology.analysis,results.ontology.analysis$effect=="Increased")

  rresults <- x.inc[order(-x.inc$counts),]
  my.title=paste(graph.title,"Increased")
  plot.inc <- .network.plot()

  x.dec <- subset(results.ontology.analysis,results.ontology.analysis$effect=="Decreased")

  rresults <- x.dec[order(x.dec$counts),]
  my.title=paste(graph.title,"Decreased")
  plot.dec <- .network.plot()

  return(list(plot.increased=plot.inc,
              plot.decreased=plot.dec,
              combined.plot=cowplot::plot_grid(plot.inc,plot.dec)))

}


plot.ont.relationships <- function(results.ontology.analysis){

  library(ggrepel)
  library(igraph)
  rresults <- results.ontology.analysis[order(-results.ontology.analysis$counts),]
  pathways <- rresults$pathway[1:15]
  pth1 <- character()
  pth2 <- character()
  pval1 <- numeric()
  pval2 <- numeric()
  enrichment <- numeric()
  w <- numeric()
  g <- character()
  counts <- numeric()
  i <- 1
  for (pathway.1 in pathways){

    count <- results[results$pathway==pathway.1,"counts"]
    pvalue <- results[results$pathway==pathway.1,"pvalues"]
    #if (count>1 & pvalue<0.05){
      genes.1 <- unlist(strsplit( results[results$pathway==pathway.1,"genes"],";"))
      for (pathway.2 in pathways){
        genes.2 <- unlist(strsplit( results[results$pathway==pathway.2,"genes"],";"))

        #if (pathway.1!=pathway.2){
          common <- intersect(genes.1,genes.2)
          if (length(common)>0){
            pth1[i] <- pathway.1
            pth2[i] <- pathway.2
            g[i] <- paste0( common,collapse = ";")
            w[i] <- length(common)
            pval1[i] <-  results[results$pathway==pathway.2,"FDR"]
            pval2[i] <-  results[results$pathway==pathway.1,"FDR"]
            enrichment[i] <- results[results$pathway==pathway.1,"enrichment"]
            counts[i] <- results[results$pathway==pathway.1,"counts"]
            i <- i+1
          }
        #}
     #}
    }

  }

  df.network <- data.frame(node.1 =pth1,
                           node.2=pth2,
                           weight=w,
                           genes=g,
                           pvalue.1 =pval1,
                           pvalue.2=pval2,
                           enrichment=enrichment,
                           counts=counts)

  df.network <- df.network[order(df.network$pvalue.1),]
  df.s <- subset(df.network,df.network$counts>2)

  xx <- reshape2::dcast(df.s, node.1~node.2, value.var = 'weight')
  xx[is.na(xx)] <- 0
  rownames(xx) <- xx$node.1

  caught.inc <- graph.incidence(xx[,2:ncol(xx)], weighted = TRUE)  #make data into a bipartite graph object
  obs.parties.all <- bipartite.projection(caught.inc)[[1]]
  obs.spp.all <- bipartite.projection(caught.inc)[[2]]

  fr.all <- layout.fruchterman.reingold(obs.spp.all)

  fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
  fr.all.df$species <- colnames(xx[,2:ncol(xx)])  ## add in the species codes

  g <- get.data.frame(obs.spp.all)
  g$from.x <- fr.all.df$V1[match(g$from, fr.all.df$species)]  #  match the from locations from the node data.frame we previously connected
  g$from.y <- fr.all.df$V2[match(g$from, fr.all.df$species)]
  g$to.x <- fr.all.df$V1[match(g$to, fr.all.df$species)]  #  match the to locations from the node data.frame we previously connected
  g$to.y <- fr.all.df$V2[match(g$to, fr.all.df$species)]

  g$pvalue.1 <- df.network$pvalue.1[match(g$from,df.network$node.1)]
  g$pvalue.2 <- df.network$pvalue.2[match(g$from,df.network$node.2)]

  fr.all.df$pvalue.1 <- df.network$pvalue.1[match(fr.all.df$species,df.network$node.1)]
  fr.all.df$pvalue.2 <- df.network$pvalue.2[match(fr.all.df$species,df.network$node.2)]

  fr.all.df$E <- df.network$enrichment[match(fr.all.df$species,df.network$node.2)]
  fr.all.df$counts <- df.network$counts[match(fr.all.df$species,df.network$node.2)]

  pplot <- ggplot() +
    geom_segment(data=g,
                 aes(x=from.x,xend = to.x, y=from.y,yend = to.y,size=(weight/10)),colour="black",linetype=1) + # add line type
    #geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes
    geom_point(data=fr.all.df,aes(x=V1,y=V2,colour=log2(E),size=counts)) +
    geom_text_repel(data=fr.all.df,aes(x=V1,y=V2,label=species,size=counts)) + # add the node labels
    #scale_colour_manual(values=c("1"="red","2"="lightblue"))+  # add colour scaling for group membership
    scale_color_gradient2(low="blue",high = "red",mid = "white")+
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
      plot.background = element_blank())

  return(pplot)

}


