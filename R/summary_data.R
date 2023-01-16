

summary.qual.data <- function(df.combi, df.design){


  df.mods <- data.frame(table( df.combi$pep_mod))
  df.mods.phospho.st <- df.mods[grepl("ST", df.mods$Var1),]

  sites <- (df.combi[,25])

  phospho.s.unique <- format(nrow((sites[grepl("S",sites$...25,fixed = T),])),big.mark = ",")
  phospho.t.unique <- format(nrow(sites[grepl("(T",sites$...25,fixed = T),]),big.mark = ",")
  phospho.y.unique <- format(nrow(sites[grepl("(Y",sites$...25,fixed = T),]),big.mark = ",")

  n.st.peptides <- sum(df.mods.phospho.st$Freq)
  df.mods.phospho.y <- df.mods[grepl("(Y)", df.mods$Var1, fixed = T),]

  n.y.peptides <- format(sum(df.mods.phospho.y$Freq),big.mark = ",")

  n.peptides.unique <- format(length(unique(unlist(sites))),big.mark = ",")

  n.peptides.total <- format(nrow(df.combi),big.mark = ",")

  nconditions <- length(unique(df.design$condition))
  nruns <- nrow(df.design)

  df.phospho <- subset(df.combi,grepl("Phospho",df.combi$pep_mod,fixed = T))
  total.phosphoproteins <- format(length(unique(df.phospho$acc_no )),big.mark = ",")
  total.proteins <- format(length(unique(df.combi$acc_no)),big.mark = ",")

  summary.text <- paste("The experiment identified", n.peptides.total,
  "peptides, of which", n.peptides.unique,
  "were unique. There were",phospho.s.unique,   "unique phosphoserine,",
  phospho.t.unique, "phosphothreonine and",
  phospho.y.unique, "phosphotyrosine sites identified and quantified in the experiment.",
  "A total of", total.proteins, "proteins were identified in the experiment of which",
  total.phosphoproteins,"were phosphorylated.",
  " \n",
  "The experiment compared",nconditions, "conditions in", nruns/nconditions,
  "replicates, requiring",nruns,"LC-MS/MS runs.")

}


summary.qual.data.proteomics <- function(df.norm, df.design,replicates=0, conditions=0){



  protein.groups <- unique(unlist(df.norm$protein.group))
  n.protein.groups <- format(length(protein.groups),big.mark = ",")
  #length(unique(unlist(df.combi[,1])))
  n.peptides.unique <- format(sum(df.norm$n.peptides),big.mark = ",")
  n.peptides.total <- format(sum(df.norm$n.psm),big.mark = ",")
  total.proteins <- format(length(unique(df.norm$acc)),big.mark = ",")

  #workout No. replicates
  if (replicates==0){
    conditions <- length(unique(df.design$condition))
    replicates <- nruns/conditions
  }
  nruns <- replicates*conditions
  data.points <- nruns*length(protein.groups)
  summary.text <- paste("The experiment produced", n.peptides.total,
                        "peptides spectral matches, of which", n.peptides.unique,
                        "were unique peptides, and belonging to", n.protein.groups, "protein groups.",
                        " \n",
                        "The experiment compared",conditions, "conditions in", replicates,
                        "replicates, requiring",nruns,"LC-MS/MS runs and producing",data.points ,"quantitative datapoints.")

}

