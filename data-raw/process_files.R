# Process species-specific data files


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

# Save the processed data in package (saves to '../data/psite_ks_species.rda')
usethis::use_data(psite_ks_species)


# Combine with previous protein and kinase-substrate sets ----

# Load list object
load("../data/protein_and_ks_sets.rda")

# Combine lists
protein_and_ks_sets <- c(protein_and_ks_sets, psite_ks_species)

# Save data in package (saves to '../data/protein_and_ks_sets.rda')
usethis::use_data(protein_and_ks_sets, overwrite = TRUE)
