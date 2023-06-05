# Load packages ----
library(devtools)
library(dplyr)


# Create species-specific PhosphoSitePlus data ----

# Load PhoshoSitePlus data
# Available from: https://www.phosphosite.org/staticDownloads
psite_Kinase_Substrate_Dataset <- read.csv("Kinase_Substrate_Dataset_psite_2023_05.csv")

# Get unique species in dataset
species_kin <- unique(psite_Kinase_Substrate_Dataset$KIN_ORGANISM)
species_sub <- unique(psite_Kinase_Substrate_Dataset$SUB_ORGANISM)
species_list <- dplyr::intersect(species_kin, species_sub)

# Create dataframe in `protools2` required format
psite_ks_species <- vector("list", length = length(species_list))
names(psite_ks_species) <- paste0("psite_", species_list)
for (species in species_list) {
  psite_ks_species[[paste0("psite_", species)]] <- psite_Kinase_Substrate_Dataset %>%
    dplyr::filter(KIN_ORGANISM == species & SUB_ORGANISM == species) %>%  # Filters matching species
    dplyr::group_by(GENE) %>%
    dplyr::summarise(
      m = dplyr::n(),
      subs = paste0(SUB_GENE, "(", SUB_MOD_RSD, ")", collapse = ";")
    ) %>%
    dplyr::mutate(
      organism = species,
      gene.db = paste0(GENE, ".pSite")
    ) %>%
    dplyr::rename(gene = GENE)
}

# Save the processed data (saves to 'data/psite_ks_species.rda')
usethis::use_data(psite_ks_species)


#####

jj <- 1

if (jj == 1) {
  edges <- "https://www.dropbox.com/s/ttmzd40mnjgh1iu/edges.csv?dl=1"
  pdts <- "https://www.dropbox.com/s/86jfnayv0qa1n2q/pdts.csv?dl=1"
  psite <- "https://www.dropbox.com/s/eb1qoofz793f4tq/psite.csv?dl=1"
  signor <- "https://www.dropbox.com/s/alpbq880emz1z2t/signor.csv?dl=1"
  reactome <- "https://www.dropbox.com/s/jdcc1355cz73mmi/reactome.csv?dl=1"

  process <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/process_ont.csv"
  myfunction <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/function_ont.csv"
  location <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/location_ont.csv"

  process_mouse <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/process_ont_mouse.csv"
  myfunction_mouse <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/function_ont_mouse.csv"
  location_mouse <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/location_ont_mouse.csv"

  process_rat <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/process_ont_rat.csv"
  myfunction_rat <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/function_ont_rat.csv"
  location_rat <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/location_ont_rat.csv"

  process_pig <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/process_ont_pig.csv"
  myfunction_pig <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/function_ont_pig.csv"
  location_pig <- "C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/location_ont_pig.csv"

  nci <- "https://www.dropbox.com/s/fe8t4nyhbljsn5y/nci.csv?dl=1"
  human.cell.markers.full <- "https://www.dropbox.com/s/a2tvvuh8kgi340v/human_cell_markers_full_lineage.csv?dl=1"
  human.cell.markers.short <- "https://www.dropbox.com/s/4njofwfs5uu3wya/human_cell_markers_short_lineage.csv?dl=1"
  mouse.cell.markers <- "https://www.dropbox.com/s/223bpbc5p24j919/mouse_cell_markers_full_lineage.csv?dl=1"
  blood.cell.markers.human1 <- "https://www.dropbox.com/s/kl2apjfnbqdw7gj/blood_cell_markers.csv?dl=1"
  blood.cell.markers.human2 <- "https://www.dropbox.com/s/o28o4fbmmuwtgs3/blood_cell_markers2.csv?dl=1"
  bone.marrow.cell.markers <- "https://www.dropbox.com/s/gl0vim65cz9804a/bone_marrow_cell_markers.csv?dl=1"
  tf.targets <- "https://www.dropbox.com/s/cuitt6gzipxaaec/TF%20target%20genes%20gtrd%20v71.csv?dl=1"

  pdt.sig <- "https://www.dropbox.com/s/wr6j8i5gosm959x/pdts_signor.csv?dl=1"

  tf.targets.omnipath <- "https://www.dropbox.com/s/hg3slk150l7zd0x/TF%20all%20targets%20omnipath.csv?dl=1"
  chromatin <- "https://www.dropbox.com/s/y2gkjn45oxy72nw/chromatin.csv?dl=1"
  selected <- "https://www.dropbox.com/s/nrraoj8hlvzap87/selected.csv?dl=1"
  ctams <- "https://www.dropbox.com/s/ev95hm4zz4c535c/ctams.csv?dl=1"
  kegg <- "https://www.dropbox.com/s/gm2821cmxarv7sx/kegg%20pathways.csv?dl=1"
  pp.markers.corr.with.cd <- "https://www.dropbox.com/s/aixtdiwhuym9vk6/PP%20Markers%20corr%20with%20CDs.csv?dl=1"

  hallmark.genes <- "https://www.dropbox.com/s/wqvnoalg2v6ufm8/hallmark%20genes%20v71.csv?dl=1"

  cd.phospho.markers <- "https://www.dropbox.com/s/aixtdiwhuym9vk6/PP%20Markers%20corr%20with%20CDs.csv?dl=1"
  pdts.reactome <- "https://www.dropbox.com/s/rvqjmzg9cq9yjpu/pdts_reactome.csv?dl=1"
  pdts.process <- "https://www.dropbox.com/s/1vqm9kky5ctzj2g/pdts_process.csv?dl=1"
  pdts.location <- "https://www.dropbox.com/s/jbroun1ojqu9rnc/pdts_location.csv?dl=1"
  pdts.nci <- "https://www.dropbox.com/s/qjvmxgp6fd5cv8x/pdts_nci.csv?dl=1"

  ks.omnipath <- "https://www.dropbox.com/s/zdnuxe8anwboks1/Kinase%20substrates%20omnipath.csv?dl=1"

  circuitries <- "https://www.dropbox.com/s/9vykdbcvy0jzh0f/Results%20circuitry%20anal.csv?dl=1"

  ctams.hijazi <- "https://www.dropbox.com/s/zksj57u60q6sg7q/CTAMS_hijazi.csv?dl=1"


  dataset.names <- c(
    "edges", # 1
    "pdts",
    "psite",
    "signor",
    "ks.omnipath", # 5
    "reactome",
    "process",
    "function",
    "location",
    "process_mouse", # 10
    "function_mouse",
    "location_mouse",
    "process_rat",
    "function_rat",
    "location_rat", # 15
    "process_pig",
    "function_pig",
    "location_pig",
    "nci", # 19
    "human.cell.markers.full", # 20
    "human.cell.markers.short",
    "mouse.cell.markers",
    "blood.cell.markers.human1",
    "blood.cell.markers.human2",
    "bone.marrow.cell.markers", # 25
    "tf.targets",
    "tf.targets.omnipath",
    "chromatin",
    "selected",
    "ctams", # 30
    "kegg",
    "cd.phospho.markers",
    "pdts.reactome",
    "pdts.process",
    "pdts.location", # 35
    "pdts.nci",
    "markers.corr.with.cd",
    "hallmark.genes",
    "pdt.sig",
    "circuitries", # 40
    "ctams.hijazi",
    "psite_mouse",  # Incorporated new species-specific PhosphoSitePlus data
    "psite_rat",
    "psite_human",
    "psite_rabbit",  # 45
    "psite_chicken",
    "psite_cow",
    "psite_pig",
    "psite_frog",
    "psite_dog",  # 50
    "psite_hamster"
  )

  datasets <- list(
    edges, # 1
    pdts, # 2
    psite,
    signor,
    ks.omnipath, # 5
    reactome,
    process,
    myfunction,
    location,
    process_mouse, # 10
    myfunction_mouse,
    location_mouse,
    process_rat,
    myfunction_rat,
    location_rat, # 15
    process_pig,
    myfunction_pig,
    location_pig,
    nci, # 19
    human.cell.markers.full, # 20
    human.cell.markers.short,
    mouse.cell.markers,
    blood.cell.markers.human1, # 23
    blood.cell.markers.human2,
    bone.marrow.cell.markers, # 25
    tf.targets,
    tf.targets.omnipath,
    chromatin, # 28
    selected,
    ctams, # 30
    kegg,
    cd.phospho.markers,
    pdts.reactome,
    pdts.process,
    pdts.location, # 35
    pdts.nci,
    pp.markers.corr.with.cd,
    hallmark.genes,
    pdt.sig,
    circuitries, # 40
    ctams.hijazi,
    psite_ks_species$psite_mouse,  # Incorporated new species-specific PhosphoSitePlus data
    psite_ks_species$psite_rat,
    psite_ks_species$psite_human,
    psite_ks_species$psite_rabbit,  # 45
    psite_ks_species$psite_chicken,
    psite_ks_species$psite_cow,
    psite_ks_species$psite_pig,
    psite_ks_species$psite_frog,
    psite_ks_species$psite_dog,  # 50
    psite_ks_species$psite_hamster
  )


  # df.datasets <- list(datasets)
  names(datasets) <- dataset.names
  dataset.list <- list()
  i <- 1
  for (d in datasets) {
    names(datasets[1])
    dataset.list[[i]] <- read.csv(datasets[[i]])
    names(dataset.list)[i] <- names(datasets[i])
    i <- i + 1
  }
  protein_and_ks_sets <- dataset.list
  k <- protein_and_ks_sets$process
  usethis::use_data_raw(name = "protein_and_ks_sets")
  usethis::use_data(protein_and_ks_sets, overwrite = TRUE)
}



uniprot.names <- read.csv("C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/uniprot_reviewed_20230102.csv")

usethis::use_data_raw(name = "uniprot.names")
usethis::use_data(uniprot.names, overwrite = TRUE)

uniprot.names.mouse <- read.csv("C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/uniprot/uniprot_mouse_reviewed_20230102.csv")

usethis::use_data_raw(name = "uniprot.names.mouse")
usethis::use_data(uniprot.names.mouse, overwrite = TRUE)

####


datasets <- protools::protein_and_ks_sets

d <- datasets$selected


j <- 10
for (d in datasets) {
  dd <- datasets[j]
  nn <- names(dd[[1]])
  ddd <- dd[[1]]
  names(dd[1])
  ddd <- subset(ddd, ddd[, 2] > 2)
  # if (!"genes" %in% nn){
  i <- 1
  for (i in 1:nrow(ddd)) {
    prots <- unique(unlist(strsplit(ddd[i, 3], ";")))
    genes <- protools::gene.names.from.accessions(prots)
    ddd$genes[i] <- paste0(genes, collapse = ";")
  }
  datasets[[j]] <- ddd
  # }
  j <- j + 1
}
protein_and_ks_sets <- datasets

usethis::use_data_raw(name = "protein_and_ks_sets")
usethis::use_data(protein_and_ks_sets, overwrite = TRUE)


df.kegg <- read.csv("C:/Users/cutill01/Dropbox/01pedro work/01_data/databases/gene sets/kegg pathways_2.csv")
i <- 1
accs <- accessions.from.gene.names(df.kegg$gene)
df.kegg$proteins <- accs

pathways <- unique(df.kegg$pathway)
nn <- length(pathways)

pathway <- character(nn)
genes <- character(nn)
proteins <- character(nn)
m <- numeric(nn)
i <- 1
for (p in pathways) {
  xx <- df.kegg[df.kegg$pathway == p, ]
  genes[i] <- paste0(xx$gene, collapse = ";")
  proteins[i] <- paste0(xx$proteins, collapse = ";")
  m[i] <- nrow(xx)
  i <- i + 1
}

df.kegg.b <- data.frame(pathways, m, proteins, genes)

df.kegg.b <- df.kegg.b[order(-df.kegg.b$m), ]

protein_and_ks_sets[["kegg"]] <- df.kegg.b

usethis::use_data_raw(name = "protein_and_ks_sets")
usethis::use_data(protein_and_ks_sets, overwrite = TRUE)

k <- protein_and_ks_sets$kegg$genes[1:10]



use_package("foreach")
use_package("doParallel")
use_package("readxl")
use_package("dplyr")
use_package("ggrepel")
use_package("stringr")
use_package("ggplot2")
use_package("ggpubr")
use_package("ggrepel")
use_package("igraph")
use_package("limma")
document()

# Reload the package: CTRL-L or
load_all()
