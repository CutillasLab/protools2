


gene.names.from.accessions <- function(accessions){
  #library(foreach)
  #library(doParallel)

  if (grepl("HUMAN", accessions)[1]==TRUE){
    prot.data <- protools::uniprot.names
  }
  if (grepl("MOUSE", accessions)[1]==TRUE){
    prot.data <- protools::uniprot.names.mouse
  }

  nr <- length(accessions)
  #print (nr)
  df.s1 <- data.frame(acc=unlist(accessions),gene="gene",stringsAsFactors = F)
  #t1 <- Sys.time()
  for (i in 1:nr ){
    acc <- NA
    #site <- as.character(phosphosite.names[i])
    #gene <- strsplit(site,"(",fixed = TRUE)[[1]][1]
    gene <- as.character(df.s1$acc[i])
    if (nchar(gene)>0){
      acc.1 <- subset(prot.data, prot.data$Entry.Name==gene)
      #if (nrow(acc.1)==1){
      acc <- as.character(acc.1$Gene.Names[1])
      #}
    }
    df.s1$gene[i] <- acc
  }
  #t2 <- Sys.time()
  #print(t2-t1)
  #print ("Gene names converted.")
  return(na.omit(df.s1$gene))
}

accessions.from.gene.names <- function(genes){


  if (grepl("HUMAN", accessions)[1]==TRUE){
    prot.data <- protools::uniprot.names
  }
  if (grepl("MOUSE", accessions)[1]==TRUE){
    prot.data <- protools::uniprot.names.mouse
  }
  nr <- length(genes)
  #print (nr)
  df.s1 <- data.frame(acc=unlist(genes),gene="gene",stringsAsFactors = F)
  #t1 <- Sys.time()
  for (i in 1:nr ){
    acc <- NA
    #site <- as.character(phosphosite.names[i])
    #gene <- strsplit(site,"(",fixed = TRUE)[[1]][1]
    gene <- as.character(df.s1$acc[i])
    if (nchar(gene)>0){
      acc.1 <- subset(prot.data,gene %in% unlist(strsplit( prot.data$Gene.Names,";")))
      #if (nrow(acc.1)==1){
      acc <- as.character(acc.1$Entry.Name[1])
      #}
    }
    df.s1$gene[i] <- acc
  }
  #t2 <- Sys.time()
  #print(t2-t1)
  #print ("Gene names converted.")
  return((df.s1$gene))
}



mycolors <- function(){
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black",
    "skyblue2",
    "palegreen2",
    "#CAB2D6", # lt purple
    "orangered", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1",  "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown", "gold1",
    "#FB9A99", "deeppink1", "royalblue" # lt pink
  )

  c23 <- c(
    "dodgerblue2",
    "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black",
    "skyblue2",
    "palegreen2",
    "#CAB2D6", # lt purple
    "orangered", # lt orange
    "gray70",
    "khaki2",
    #"maroon",
    #"orchid1",
    "blue1", "steelblue4",
    "darkturquoise",
    "green1", "yellow4", "yellow3",
    "darkorange4", "brown", "gold1",
    "#FB9A99", "deeppink1", "royalblue" # lt pink
  )
  return(list(C23=c23,C25=c25))
}
