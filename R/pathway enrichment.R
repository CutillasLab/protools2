
available_protein_and_ks_datasets <- function(){
  print(names(protools::protein_and_ks_sets))
}


expand.phosphopeptide.dataset <- function(df.fold){

  pps <- unique(unlist(strsplit( df.fold[,1],";")))
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  t1 <- Sys.time()
  dfx <- foreach (p = pps,.combine = "rbind")%dopar% {
    data.frame(p, df.fold[grepl(p, df.fold[,1],fixed = T),2:ncol(df.fold)])
  }
  stopCluster(cl)
  t2 <- Sys.time()
  return(dfx)
}

ksea <- function(df.fold,ks_db=c("pdts","psite"), graph.heading="",pval.cut.off=0.05){

  library(ggrepel)

  yy <- expand.phosphopeptide.dataset(df.fold)
  cores=detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  xx <- foreach(db = ks_db,.combine = "rbind")%do%{
    x <- kinase.substrate.enrichment(dfx= yy,ks_db =  db)
    x
  }
  stopCluster(cl)

  df.k.up <- subset(xx, xx$pvalues<pval.cut.off & xx$zscores>0)
  df.k.do <- subset(xx, xx$pvalues<pval.cut.off & xx$zscores<0)

  plot.v.ksea.p <- ggplot(xx,aes(x=zscores,y=-log10(pvalues)))+
    geom_point()+
    geom_point(data = df.k.up ,aes(x=zscores,y=-log10(pvalues),size=m),color="red")+
    geom_label_repel(data = df.k.up ,aes(x=zscores,y=-log10(pvalues),label=kinase_dataset),color="red")+
    geom_point(data = df.k.do ,aes(x=zscores,y=-log10(pvalues),size=m),color="blue")+
    geom_label_repel(data = df.k.do ,aes(x=zscores,y=-log10(pvalues),label=kinase_dataset),color="blue")+
    labs(title = graph.heading, subtitle = "Kinase Substrate Enrichment Analysis (KSEA)", caption = "Kinase activities quantified from changes in substrate phosphorylation")+
    geom_vline(xintercept = 0, linetype=2)+
    theme_bw()

  df.k.up <- subset(xx, xx$qvalue<pval.cut.off & xx$zscores>0)
  df.k.do <- subset(xx, xx$qvalue<pval.cut.off & xx$zscores<0)

  plot.v.ksea.q <- ggplot(xx,aes(x=zscores,y=-log10(qvalue)))+
    geom_point()+
    geom_point(data = df.k.up ,aes(x=zscores,y=-log10(qvalue),size=m),color="red")+
    geom_label_repel(data = df.k.up ,aes(x=zscores,y=-log10(qvalue),label=kinase_dataset),color="red")+
    geom_point(data = df.k.do ,aes(x=zscores,y=-log10(qvalue),size=m),color="blue")+
    geom_label_repel(data = df.k.do ,aes(x=zscores,y=-log10(qvalue),label=kinase_dataset),color="blue")+
    labs(title = graph.heading, subtitle = "Kinase Substrate Enrichment Analysis (KSEA)", caption = "Kinase activities quantified from changes in substrate phosphorylation")+
    geom_vline(xintercept = 0, linetype=2)+
    theme_bw()

  #plot.v.ksea
  return(list(ksea.data=xx, ksea.volcano.plot.pvals=plot.v.ksea.p, ksea.volcano.plot.qvals=plot.v.ksea.q))

 }

kinase.substrate.enrichment <- function(dfx, ks_db,is.ksea=TRUE){



  # df.fold == dataset of fold changes or
  # ks_db == database of kinase-substrate relationships
  #       possibilities are "edges", "ctams" "pdts", "pSite"
  #
  # returns 6 data frames:
  #         results.zscores,
  #         results.distance
  #         results.pvalues,
  #         results.m,
  #         results.q
  #         results.sites
  dfx <- data.frame(dfx)
  rownames(dfx) <- dfx[,1]
  if (is.ksea==TRUE){
    column.with.priors <- 3
    }else{
    #Convert phosphopeptides to genes for phosphoproteomics data when matching to proteomics sets
    x1 <-  grepl("(", as.character(dfx[,1]),fixed=T)
    if (TRUE %in% x1){ # True if phosphoproteomics data, false if proteomics data
      dfx[,1] <- unlist(lapply(dfx[,1],
                                        function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1]))
      column.with.priors <- "genes"
    }
  }

  dfx[,1] <- gsub("..",");",dfx[,1],fixed = T)
  dfx[,1] <- gsub(".","(",dfx[,1],fixed = T)
  nc <- ncol(dfx)
  df.ks <- protools::protein_and_ks_sets[[ks_db]]
  nr <- nrow(df.ks)
  results.pvalues <-numeric(nr)
  results.zscores <-numeric(nr)
  results.q <-numeric(nr)
  results.m <-numeric(nr)
  results.distance <-numeric(nr)
  results.sites <-character(nr)
  kinases <- character(nr)
  r=1
  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,column.with.priors])
        ss <- c(unlist(strsplit(substrates,";")))
        start.time <- Sys.time()

        df.xx <- subset(dfx,dfx[,1] %in% ss)# paste(ss,";",sep=""))

        sites <-paste(unlist(rownames(df.xx)),collapse=";")


        if (nrow(df.xx)>2){
          df.xx$prot.group <- kinase
          sites.x <- data.frame(site=df.xx[,1])
          sites.x$kinase <- kinase
          #all.sites <- rbind(all.sites,sites.x)
          c=2
          ds <- numeric(nc-1)
          pvals <- numeric(nc-1)
          zscores <- numeric(nc-1)
          ms<- numeric(nc-1)
          qs <- numeric(nc-1)

            values.all <- na.omit(as.numeric(subset(dfx[,c], dfx[,c]!=0)))
            myvalues <- na.omit(as.numeric(subset(df.xx[,c],df.xx[,c]!=0)))
            pval <- 1
            tryCatch({
              myks <- ks.test(values.all,myvalues)
              pval <- myks$p.value
            }, error=function(e){}
            )

            m <- 0
            q <- 0
            m <- nrow(df.xx)
            mysd <- sd(values.all,na.rm = T)
            mymedian <- median(myvalues, na.rm = T)
            mymedian.all <- median(values.all,na.rm = T)
            sd.all <- sd(values.all)

          results.distance[r] <- mymedian-mymedian.all
          results.pvalues[r] <- pval
          results.zscores[r] <- ((mymedian-mymedian.all)*((sqrt(m)))/mysd)
          results.m[r] <- m
         # results.q[r] <- qs
          results.sites[r] <- as.character(sites)
          kinases[r] <- kinase

        }
      }
    }
  }

  xx <- data.frame(kinases,
                   zscores=results.zscores,
                   pvalues=results.pvalues,
                   m=results.m,
                   #q=results.q,
                   distance=results.distance,
                   sites=results.sites,
                   kinase_dataset =paste0(kinases,"_",ks_db))

  xx <- subset(xx,xx$m>1)

  #xx <- xx[order(xx$zcores),]
  xx$qvalue <- p.adjust(xx$pvalues,method = "fdr")
  #head(xx[order(xx$pvalues),],n=20)
  return(xx)
  ##################################
}


pathway_enrichment_matrix <- function(list_of_comparisons_from_limma,
                                      qval_cutoff=0.25,
                                      fold_cutoff=0,
                                      pathway_pval_cutoff=0.05,
                                      pathway_diff_cutoff=0,
                                      prot_dbs=c("kegg","hallmark.genes","nci","process"),
                                      n_pathways=50){

  #list_of_comparisons_from_limma <- comparisons

  ##list.of.vps <- list()
  #list.of.delta.enrichment.plots <- list()
  #list.of.pathway.enrichment.plots <- list()

  delta.vals <- list()
  top.pths <- character()
  names(list_of_comparisons_from_limma)
  i <- 1
  for (comp in list_of_comparisons_from_limma){
    dd <- list_of_comparisons_from_limma[i]
    names(dd)
    # Identify differences in comparisons
    suppressMessages(
      ddd <- identify_differences_in_comparison_plus_volcano(dd[[1]],
                                                             fold.cutoff = fold_cutoff,
                                                             qval.cutoff = qval_cutoff,
                                                             graph.header=names(dd))
    )

    decreased.peptides <- ddd$df.decreased$sites
    increased.peptides <- ddd$df.increased$sites

    ppx <- data.frame(pathway=NA,
               delta.enrichment=0,
               pvalue=1,
               delta.count=0,
               sample=names(dd))

    if (length(decreased.peptides)>0 & length(increased.peptides)>0){
      pe <- protools::pathway_enrichment(increased.peptides = increased.peptides,
                                         decreased.peptides = decreased.peptides,
                                         background.data =  dd[[1]],
                                         graph.heading = names(dd),
                                        prot_dbs=prot_dbs)

      if (length(pe)>1){
        df.delta <- pe$delta_pathway_enrichment_data
        if (nrow(df.delta)>0){
          df.delta$best.qvalue <- min(df.delta$FDR.x,df.delta$FDR.y)
          ppx <- data.frame(pathway=df.delta$pathway,
                            delta.enrichment=df.delta$delta.enrichment.up.vs.down,
                            pvalue=df.delta$best.qvalue,
                            delta.count=df.delta$delta.counts,
                            sample=names(dd)
                            )
          ppx <- subset(ppx,ppx$counts>n_counts_cutoff & ppx$pvalue<pathway_pval_cutoff)

          if (nrow(ppx)>0){
            ppx <- subset(ppx, ppx$qvalue<qval_cutoff & abs(ppx$delta.enrichment.up.vs.down)>pathway_diff_cutoff)
            ppx$alpha <- ppx$delta.enrichment*(-log10(ppx$pvalue))

            ppx <- ppx[order(ppx$alpha),]
            top.pths <-c(top.pths, unlist(ppx$pathway[1:n_pathways]))
            ppx <- ppx[order(-ppx$alpha),]
            top.pths <-c(top.pths, unlist(ppx$pathway[1:n_pathways]))
          }
        }
      }
    }

    delta.vals[[i]]   <- ppx

    i <- i+1
  }



  df.delta.counts <- do.call(rbind.data.frame, delta.vals)

  top.pths <- na.omit(top.pths)

  sig.pths <-(top.pths)
  xsig.pths <- data.frame(table(sig.pths))
  xsig.pths <- xsig.pths[order(-xsig.pths$Freq),]

  ss <- na.omit(df.delta.counts[df.delta.counts$pathway %in% xsig.pths$sig.pths[1:n_pathways],])

  a <- ss$pathway
  a <- ifelse(nchar(a) > 35, paste0(strtrim(a, 35), '...'), a)
  ss$pathway <- a

  ssm <- reshape2::dcast(ss,pathway~sample, value.var = "delta.enrichment")
  avs <- data.frame(pathway = ssm$pathway,  av=apply(ssm[,2:ncol(ssm)],1,mean))
  avs <- avs[order(avs$av),]



 pp <-  ggplot(na.omit(ss),aes(x=sample,y=pathway))+
    geom_point(aes(color=delta.enrichment,size=-log10(pvalue)))+
    labs(title = "Pathway Enrichement Analysis",
         subtitle = "Methods = Fold enrichment and hypergeometric test",
         color="delta(E)", size="-log10(q)", x= "Sample",y="")+
   scale_color_gradient2(low="blue", mid="white", high = "red")+
   scale_y_discrete(limits=avs$pathway)+
    theme_bw()
  #head(df.all)

 return(list(pathway_enrichment_data=df.delta.counts,
             pathway_enrichment_data_significant_values=ss,
             plot_enrichment=pp))

}



pathway_enrichment_matrix_by_page_and_ks <- function(list_of_comparisons_from_limma,
                                                     qval_cutoff=1,
                                                     fold_cutoff=0,
                                                     prot_dbs=c("kegg","hallmark.genes","nci","process"),
                                                     n_pathways=50,
                                                     n_counts_cutoff =1,
                                                     pval_cutoff= 0.1){


  library(foreach)
  library(doParallel)
  library(ggplot2)

  zscores <- list()

  top.pths <- character()
  #names(list_of_comparisons_from_limma)
  i <- 1
  for (comp in list_of_comparisons_from_limma){

    dd <- list_of_comparisons_from_limma[i]
    cores <- detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    xx <- foreach(db = prot_dbs,.combine = "rbind")%dopar%{
      x <- protools::kinase.substrate.enrichment(dfx= dd[[1]],ks_db =  db, is.ksea=FALSE)
      x
    }
    stopCluster(cl)

    ppx <- data.frame(pathway=xx$kinases, zscore=xx$zscores,pvalue=xx$pvalues,
                               qvalue=xx$qvalue,counts=xx$m, alpha=xx$zscores*(-log10(xx$pvalues)),
                               sites=xx$sites,  sample=names(dd))



    zscores[[i]] <- ppx

    ppx <- subset(ppx,ppx$counts>n_counts_cutoff & ppx$pvalue<pval_cutoff)
    ppx <- subset(ppx, ppx$qvalue<qval_cutoff & abs(ppx$zscore)>fold_cutoff)
    ppx <- ppx[order(ppx$alpha),]
    top.pths <-c(top.pths, unlist(ppx$pathway[1:25]))
    ppx <- ppx[order(-ppx$alpha),]
    top.pths <-c(top.pths, unlist(ppx$pathway[1:25]))
    i <- i+1
  }

  df.zscores <- do.call(rbind.data.frame, zscores)
  top.pths <- na.omit(top.pths)
  sig.pths <-(top.pths)
  xsig.pths <- data.frame(table(sig.pths))
  xsig.pths <- xsig.pths[order(-xsig.pths$Freq),]

  ss <- na.omit(df.zscores[df.zscores$pathway %in% xsig.pths$sig.pths[1:n_pathways],])
  a <- ss$pathway
  a <- ifelse(nchar(a) > 35, paste0(strtrim(a, 35), '...'), a)
  ss$pathway <- a
  ssm <- reshape2::dcast(ss,pathway~sample, value.var = "zscore")
  avs <- data.frame(pathway = ssm$pathway,  av=apply(ssm[,2:ncol(ssm)],1,mean, na.rm=T))
  avs <- avs[order(avs$av),]

  pp <-  ggplot(na.omit(ss),aes(x=sample,y=pathway))+
    geom_point(aes(color=zscore,size=-log10(pvalue)))+
    labs(title = "Pathway Enrichement Analysis",
         subtitle = "Methods = PAGE (zscore) vs KS (pvalue)",
         color="zscore", size="-log10(P)", x= "Sample",y="")+
    scale_color_gradient2(low="blue", mid="white", high = "red")+
    scale_y_discrete(limits=avs$pathway)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  #pp

  #head(df.zscores)

  return(list(pathway_enrichment_data=df.zscores,
              pathway_enrichment_data_significant_values=ss,
              plot_enrichment=pp))

}

pathway_enrichment <- function(increased.peptides,
                               decreased.peptides,
                               background.data,
                               prot_dbs=c("kegg","hallmark.genes","nci","process"),
                               graph.heading=""
                               ){

  library(foreach)
  library(doParallel)
  background.list <- background.data$protein
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  t1 <- Sys.time()
  enrich.up <- foreach(db = prot_dbs, .combine = "rbind")%dopar%{

    e <- protools::enrichment.from.list(list.of.peptides=increased.peptides,
                                        background.list,
                                        prot_db = db,
                                        is.ksea = F)
    if (nrow(e)>0){
      e$prot_db <- db
      e
    }
  }
  stopCluster(cl)
  t2 <- Sys.time()
  #print(t2-t1)
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  t1 <- Sys.time()
  enrich.do <- foreach(pp = prot_dbs, .combine = "rbind")%dopar%{

    e <- protools::enrichment.from.list(list.of.peptides =  decreased.peptides,background.list, prot_db=pp, is.ksea = F)
    if (nrow(e)>0){
      e$prot_db <- pp
      e
    }
  }
  stopCluster(cl)
  t2 <- Sys.time()
  #print(t2-t1)

  if (length(enrich.do)>0 & length(enrich.up)>0){

    enrich.up <- enrich.up[ order(enrich.up$pvalues),]

    a <- enrich.up$pathway
    a <- ifelse(nchar(a) > 50, paste0(strtrim(a, 50), '...'), a)
    enrich.up$pathway <- a

    enrich.do <- enrich.do[ order(enrich.do$pvalues),]

    a <- enrich.do$pathway
    a <- ifelse(nchar(a) > 50, paste0(strtrim(a, 50), '...'), a)
    enrich.do$pathway <- a

    # enrich.up <- subset(enrich.up,enrich.up$counts>1)

    enrich.do$effect <- "Decreased"
    enrich.up$effect <- "Increased"
    enrich.up.do <- rbind.data.frame(enrich.do[1:25,],enrich.up[1:25,])

    ##### Differential enrichment up vs down ########################

    diff.enrich <- merge.data.frame(enrich.do,enrich.up,by="pathway")

    diff.enrich$delta.enrichment.up.vs.down <- log2(diff.enrich$enrichment.y/diff.enrich$enrichment.x)
    diff.enrich$delta.pval.up.vs.down <- -log10(diff.enrich$pvalues.y)-(-log10(diff.enrich$pvalues.x))
    diff.enrich$delta.counts <- diff.enrich$counts.y-diff.enrich$counts.x

    diff.enrich <- diff.enrich[order(-diff.enrich$delta.enrichment.up.vs.down),]

    dif.e.up <- diff.enrich[1:15,]

    diff.enrich <- diff.enrich[order(diff.enrich$delta.enrichment.up.vs.down),]
    dif.e.do <- diff.enrich[1:15,]

    ######### Enrichment both ###########################################


    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    t1 <- Sys.time()
    enrich.combined <- foreach(db = prot_dbs, .combine = "rbind")%dopar%{

      e <- protools::enrichment.from.list(list.of.peptides=c(increased.peptides,decreased.peptides),
                                          background.list,
                                          prot_db = db,
                                          is.ksea = F)
      if (nrow(e)>0){
        e$prot_db <- db
        e
      }
    }
    stopCluster(cl)
    t2 <- Sys.time()
    #print(t2-t1)

    enrich.combined <- enrich.combined[ order(enrich.combined$pvalues),]

    a <- enrich.combined$pathway
    a <- ifelse(nchar(a) > 50, paste0(strtrim(a, 50), '...'), a)
    enrich.combined$pathway <- a

    enrich.combined <- enrich.combined[ order(enrich.combined$pvalues),]

    ##### Plot results #################################################
    plot.diff <- ggplot(rbind.data.frame(dif.e.do,dif.e.up),aes(x=delta.enrichment.up.vs.down,y=reorder(pathway,
                                                                                                        delta.enrichment.up.vs.down)))+
      geom_point(aes(size=abs(delta.counts),color=delta.pval.up.vs.down))+
      labs(title = "Pathway Enrichement Analysis", subtitle = paste0("Changed in ", graph.heading),y="",
           size="delta \ncounts", color="delta \npvalue")+
      scale_color_gradient2(low="blue",high = "red",mid = "yellow")+
      geom_vline(xintercept = 0, linetype=2)+
      theme_bw()

    plot.pathways.enrichment.by.counts <- ggplot(na.omit(enrich.up.do),aes(x=counts,y=reorder(pathway,counts)))+
      geom_point(aes(size=log2(enrichment),color=-log10(FDR)))+
      labs(title = "Pathway Enrichement Analysis", subtitle = paste0("Changed in ", graph.heading),
           size="log2(E)", color="-log10(q)", x= "# Proteins",y="")+
      scale_color_gradient(low="orange",high = "red")+
      facet_wrap(~effect,scales = "free")+
      theme_bw()


    plot.pathways.enrichment.by.counts.combined <- ggplot(na.omit(enrich.combined),
                                                          aes(x=counts,y=reorder(pathway,counts)))+
      geom_point(aes(size=log2(enrichment),color=-log10(FDR)))+
      labs(title = "Pathway Enrichement Analysis", subtitle = "Combined analysis using increased and decreased peptides",
           size="log2(E)", color="-log10(q)", x= "# Proteins",y="")+
      scale_color_gradient(low="orange",high = "red")+
      facet_wrap(~effect,scales = "free")+
      theme_bw()



    plot.pathways.enrichment.by.enrichment <- ggplot(na.omit(enrich.up.do),aes(x=log2(enrichment),
                                                                               y=reorder(pathway,enrichment)))+
      geom_point(aes(size=counts,color=-log10(FDR)))+
      labs(title = "Pathway Enrichement Analysis", color="-log10(q)",
           subtitle = paste0("Changed in ", graph.heading),y="")+
      scale_color_gradient(low="orange",high = "red")+
      facet_wrap(~effect,scales = "free")+
      theme_bw()

    #write.csv(rbind.data.frame(enrich.do[1:50,],enrich.up[1:50,]),
    #          paste0("Results pathway analysis for ", graph.heading,".csv"))

    ##### look at proteins in significant ontologies #############

    #enrich.up[order(-enrich.up$enrichment),]

    uniprot.data <- protools::uniprot.names

    enrich.data <- enrich.up

    .look.at.proteins.in.pathways <- function(enrich.data){

    mypathways <- list()

    j <- 1
    mypathways <- foreach(j = 1:10,.combine = "rbind")%do%{
      ontology <- enrich.data[j,]$pathway
      prots <- unlist(strsplit( enrich.data[j,"proteins"],";"))
      #subset(background.data,background.data$protein %in% paste0(prots,";"))
      #x <- background.data$protein
      prot <- prots[1]

      dfxx <- foreach(prot =prots,.combine = "rbind")%do%{
        vv <- background.data[grepl(prot,background.data$protein,fixed = T),]
        vv <- vv[order(vv$pvalues),]
        vv[1,]
      }

      prots <- unique(unlist(strsplit( dfxx$protein,";")))

      df.names <- foreach(ppp =unlist(dfxx$protein),.combine = "rbind")%do%{
        prot <- unlist(strsplit( ppp,";"))
        fff <- foreach(pp = prot,.combine="rbind")%do%{
          p <- strsplit(pp,"(",fixed = T)[[1]][1]
          nn <- uniprot.data[grepl(p, uniprot.data$gene.name,fixed = T),]
          accs <- paste0(nn$acc,collapse = ";")
          names <- paste0(nn$protein,collapse = ";")
          data.frame(protein=ppp,accs,names)
        }
        data.frame(protein=ppp, acc=paste0(fff$accs,collapse = ";"), name=paste0(fff$names,collapse = ";"))
      }

      df.comb <- merge.data.frame(dfxx,df.names,by="protein")

      df.comb$ontology <- ontology
      df.comb
    }

    colnames(mypathways)[2] <- "fold"

    mypathways$name.short <- paste(stringr::str_trunc(mypathways$protein, 15),
                                stringr::str_trunc(mypathways$name, 40))

    p <- strsplit(pp,"(",fixed = T)[[1]][1]

    ff <- data.frame(table(mypathways$ontology))
    ff <- ff[order(ff$Freq),]
    plot.all <- ggplot(mypathways,aes(x=ontology,y=name.short))+
      scale_color_manual(values=protools::mycolors()$C23 )+
     # scale_x_discrete(limits=rev(ff$Var1))+
      geom_point(aes(size=fold,color=ontology))  +
      theme_bw()+
      theme(axis.text.x = element_blank())

    list.ont <- list()
    ii <- 1
    for (ont in unique(mypathways$ontology)){

      h <- mypathways[mypathways$ontology==ont,]
      p <- ggplot(h,aes(y=reorder(name.short,fold),x=fold))+
        geom_point(aes(color=-log(pvalues),size=-log(pvalues)))+
        #ggtitle(ont)+
        labs(title = ont,color="-log10(p)",size="",y="")+
        scale_color_gradient2(low="seagreen",high = "purple4")+
        theme_bw()+
        theme(plot.title = element_text(hjust=1))
      list.ont[[ii]] <-  list(data.table=h, dot.pot=p)
      names(list.ont)[[ii]] <- ont
      ii <- ii+1
      }

    return(list(combined.plot=plot.all,list.of.plots.and.data=list.ont))

    }


    prots.in.increased.pathways <- .look.at.proteins.in.pathways(enrich.up)
    prots.in.decreased.pathways <- .look.at.proteins.in.pathways(enrich.do)

    }else{
    return(0)
  }


  return(list(pathway_enrichment_data=rbind.data.frame(enrich.do[1:50,],enrich.up[1:50,]),
              delta_pathway_enrichment_data=rbind.data.frame(diff.enrich),
              plot_delta_enrichment=plot.diff,
              plot.pathways.enrichment.by.counts=plot.pathways.enrichment.by.counts,
              plot.pathways.enrichment.by.enrichment=plot.pathways.enrichment.by.enrichment,
              plot.pathways.enrichment.by.counts.combined=plot.pathways.enrichment.by.counts.combined,
              prots.in.increased.pathways=prots.in.increased.pathways,
              prots.in.decreased.pathways=prots.in.decreased.pathways)
  )

}

enrichment.from.list <- function(
        list.of.peptides,
         background.list,
         prot_db= c("kegg","hallmark.genes"), is.ksea=F){

  # Returns enrichment of phosphosites, proteins or transcripts in pathways ontologies or kinase-phosphosite relationsips
  #
  # list of peptides is the list of proteins, phosphosites or transcripts to be searched (transcripts needt o be converted to Uniprot names)
  # Accepts Uniprot names or phosphosite names in the form gene(xy), were gene = uniprot gene name; x=S, T or Y; y= amino acid position of phosphorylation
  # background.list is the background list of of proteins, phosphosites or transcripts
  # prot_db = any suitable ontology or kinase-substrate database
  # is.ksea = change to TRUE for analysing phosphoproteomics data against kinase-substrates
  # draw.tables.and.plots = if TRUE draws volcano plots of the enrichment data, requires ggplot2
  #
  # Enrichment = (q/k)/(j/m), where:
  #   q = peptides/proteins in list.of.peptides with a match in the ontology/K-S relationship dataset
  #   k = number of peptides/proteins in list.of.peptides
  #   j = peptides/proteins in background.list with a match in the ontology/K-S relationship dataset
  #   m = number of peptides/proteins in background.list
  #
  # p-values are calculated using the hypergeometric function and adjusted by FDR method
  #

  ll <- list.of.peptides[!grepl("_no_mod",list.of.peptides,fixed = T)]
  ll <- ll[!grepl("_Phospho (",ll,fixed=T)]
  list.of.peptides <- strsplit(ll,";")
  original.list.of.peptides <- list.of.peptides
  background.list <- strsplit(background.list,";")
  column.with.priors <- "genes"

  df.peptides <- data.frame(peptides=unlist(original.list.of.peptides),
                            proteins=unlist(list.of.peptides))

  if (is.ksea==TRUE){
    column.with.priors <- 3
  }else{

    #Convert phosphopeptides to genes for phosphoproteomics data when matching to proteomics sets
    x1 <-  grepl(")", as.character(list.of.peptides),fixed=T)
    if (TRUE %in% x1){ # True if phosphoproteomics data, false if proteomics data
      list.of.peptides <- unlist(lapply(df.peptides$proteins,
             function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1]))
      df.peptides$proteins <- list.of.peptides

      background.list <- unlist(lapply(background.list,
                                        function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1]))
    }else{
      column.with.priors <- "genes"
    }
  }

  ####################

  df.ks <- protools::protein_and_ks_sets[[prot_db]]
  df.ks <- df.ks[order(-df.ks[,2]),]
  nr <- nrow(df.ks)

  pathway <- character(nr)
  pvalues <- numeric(nr)
  enrichment <- numeric(nr)
  counts <- numeric(nr)
  counts.bg <- numeric(nr)
  FDR <- numeric(nr)
  proteins <- character(nr)
  genes <- character(nr)

  m <- length(background.list)
  n <- length(list.of.peptides)

  bg.size <- numeric(nr)
  data.size <- numeric(nr)
  r=1
  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,column.with.priors])
        ss <- c(unlist(strsplit(substrates,";")))
        if (is.ksea==TRUE | is.ksea==T){
          ss <- paste(ss,";",sep = "")
        }
        start.time <- Sys.time()
        k <- n#length(ss)
        q <- length(intersect(ss,list.of.peptides))
        j <- length(intersect(ss,background.list))
        if (q>1 & j>1){
          prots <- intersect(ss,list.of.peptides)
          peptides <- df.peptides[df.peptides$proteins %in% prots,"peptides"]
          pvalue <- 1-phyper(q-1,j,m-j,k,lower.tail = TRUE, log.p = F)
          pathway[r] <- as.character(kinase)
          pvalues[r] <- pvalue
          enrichment[r] <- round((q/k)/(j/m), digits = 2)
          counts[r] <- q
          data.size[r] <- k
          counts.bg[r] <- j
          bg.size[r] <- m
          genes[r] <- paste(peptides,collapse = ";")
          proteins[r] <- paste(prots,collapse = ";")
        }
      }
      r=r+1
    }
  }
  results <- data.frame(pathway, pvalues,enrichment,counts, data.size,counts.bg,bg.size, proteins,genes)
  results <- na.omit(results)
  if (nrow(results)>0){

    results <- subset(results, pvalues!=0 & counts>0)
    results <- results[order(-results$enrichment),]
    results$FDR <- p.adjust(results$pvalues, method = "BH")
  }
  return(results)
}
