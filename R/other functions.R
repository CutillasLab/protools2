


gene.names.from.accessions <- function(accessions, species="HUMAN"){

  library(foreach)
  library(doParallel)


  if (species == "HUMAN"){
    prot.data <- protools2::uniprot.names
  }
  if (species == "MOUSE"){
    prot.data <- protools2::uniprot.names.mouse
  }
  prot.data$Entry.Name <- paste0(";",prot.data$Entry.Name,";")
  nr <- length(accessions)
  #print (nr)
  df.s1 <- data.frame(acc=paste0(";", unlist(accessions),";"),gene="gene",stringsAsFactors = F)
  #t1 <- Sys.time()

  cl <- makeCluster(detectCores(logical = TRUE))
  registerDoParallel(cl)
  genes <- foreach(i = 1:nr, .combine = rbind)%dopar%{
    acc <- NA
    gene <- as.character(df.s1$acc[i])
    if (nchar(gene)>0){
      acc.1 <- prot.data[grepl(gene,prot.data$Entry.Name,fixed = T),]
      gene <- as.character(acc.1$Gene.Names[1])
    }
    gene
  }
  stopCluster(cl)

  return(genes)
}

accessions.from.gene.names <- function(genes, species="HUMAN"){


  library(foreach)
  library(doParallel)


  if (species == "HUMAN"){
    prot.data <- protools2::uniprot.names
  }
  if (species == "MOUSE"){
    prot.data <- protools2::uniprot.names.mouse
  }
  prot.data$genes <- paste0(";",prot.data$Gene.Names,";")
  nr <- length(genes)
  #print (nr)
  df.s1 <- data.frame(acc=paste0(";", unlist(genes),";"),gene="gene",stringsAsFactors = F)
  #t1 <- Sys.time()

  cl <- makeCluster(detectCores(logical = TRUE))
  registerDoParallel(cl)
  accs <- foreach(i = 1:nr, .combine = rbind)%dopar%{
  #for (i in 1:nr ){
    acc <- NA
    gene <- as.character(df.s1$acc[i])
    if (nchar(gene)>0){
      acc.1 <- prot.data[grepl(gene,prot.data$genes,fixed = T),]
      acc <- as.character(acc.1$Entry.Name[1])
      #}
    }
    acc
  }
  stopCluster(cl)

  return(accs)
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
