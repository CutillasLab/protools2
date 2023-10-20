# Process species-specific data files

# NOTE: Script assumes working directory is protools2 project root,
# i.e. the parent directory of data-raw (the location of this process_files.R script)

# Load packages ----
library(devtools)
library(dplyr)


# Create species-specific PhosphoSitePlus data ----

# Load PhoshoSitePlus data
# Available from: https://www.phosphosite.org/staticDownloads
psite_Kinase_Substrate_Dataset <- read.csv("data-raw/Kinase_Substrate_Dataset_psite_2023_05.csv")

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

# Save the processed data in package (saves to '../data/psite_ks_species.rda')
usethis::use_data(psite_ks_species, overwrite = TRUE)


# Combine with previous protein and kinase-substrate sets ----

# Load list object
load("../data/protein_and_ks_sets.rda")

# Combine lists
protein_and_ks_sets <- c(protein_and_ks_sets, psite_ks_species)

# Save data in package (saves to '../data/protein_and_ks_sets.rda')
usethis::use_data(protein_and_ks_sets, overwrite = TRUE)



# Uniprot ----

# Create Uniprot human (Homo sapiens) data ----
# (taxonomy_id:9606 AND reviewed:true)
uniprot_human <- readr::read_tsv("data-raw/uniprot_human.tsv")

# Make column names safe
colnames(uniprot_human) <- make.names(colnames(uniprot_human))

# Save data in package (saves to '../data/uniprot_human.rda')
usethis::use_data(uniprot_human, overwrite = TRUE)


# Create Uniprot mouse (Mus musculus) data ----
# (taxonomy_id:10090 AND reviewed:true)
uniprot_mouse <- readr::read_tsv("data-raw/uniprot_mouse.tsv")

# Make column names safe
colnames(uniprot_mouse) <- make.names(colnames(uniprot_mouse))

# Save data in package (saves to '../data/uniprot_mouse.rda')
usethis::use_data(uniprot_mouse, overwrite = TRUE)


# Create Uniprot rat (Rattus norvegicus) data ----
# (taxonomy_id:10116 AND reviewed:true)
uniprot_rat <- readr::read_tsv("data-raw/uniprot_rat.tsv")

# Make column names safe
colnames(uniprot_rat) <- make.names(colnames(uniprot_rat))

# Save data in package (saves to '../data/uniprot_rat.rda')
usethis::use_data(uniprot_rat, overwrite = TRUE)


# Create Uniprot pig (Sus scrofa) data ----
# (taxonomy_id:9823 AND reviewed:true)
uniprot_pig <- readr::read_tsv("data-raw/uniprot_pig.tsv")

# Make column names safe
colnames(uniprot_pig) <- make.names(colnames(uniprot_pig))

# Save data in package (saves to '../data/uniprot_pig.rda')
usethis::use_data(uniprot_pig, overwrite = TRUE)



# Manage package build ----
package_imports <- c(
  "knitr",
  "devtools",
  "openxlsx",
  "readxl",
  "kableExtra",
  "flextable",
  "gplots",
  "grid",
  "gridExtra",
  "cowplot",
  "reshape2",
  "viridisLite",
  "webshot",
  "webshot2",
  "plotly",
  "heatmaply",
  "tools",
  "ggrepel",
  "ggpubr",
  "igraph",
  "data.table",
  "doParallel",
  "foreach",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "tibble",
  "stringr",
  "forcats",
  "dtplyr",
  "limma"
)

lapply(
  package_imports,
  function(pkg) {
    usethis::use_package(package = pkg, type = "Imports")
  }
)
