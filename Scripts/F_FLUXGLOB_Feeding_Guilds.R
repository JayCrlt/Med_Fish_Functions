#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("rphylopic") ; library("FNN") ; library("pbapply")

## Charge from previous scripts
load("Outputs/FLUXGLOB/dat_proc/All_Species_FLUXGLOB.RData")
load("Data/Datos Atlántico/community_and_traits.RData")
## Functions
source("Scripts/00_functions_script.R")

## Generate Feeding habits 
## (https://www.fao.org/fishery/docs/CDrom/FAO_Training/FAO_Training/General/x6709e/.!17816!x6709e10.htm)
# Three sources to merge
# 1) Carlot et al. 2027
MEDFLUX <- readxl::read_xlsx("Data/Guilds_MED.xlsx") |> mutate(SPECIES = str_replace_all(SPECIES, "\u00A0", " ")) |> 
  rename(Species_FLUXGLOB = SPECIES, Feeding_mode = Trophic_category)

# 2) B-Useful
dTraits = dTraits |> dplyr::select(taxon, feeding.mode) |> rename(Species_FLUXGLOB = taxon, Feeding_mode = feeding.mode)

# 3) RFishBase
Feeding_habits = rfishbase::ecology() |> dplyr::select(SpecCode, FeedingType) |> 
  dplyr::left_join(tibble::tribble(~FeedingType, ~Feeding_mode, "browsing on substrate", "Benthivorous",
    "feeding on a host (parasite)", NA_character_, "feeding on dead animals (scavenger)", "Benthivorous",
    "feeding on the prey of a host (commensal)", NA_character_, "filtering plankton", "Planktivorous",
    "grazing on aquatic plants", "Herbivorous", "hunting macrofauna (predator)", NA_character_, "other",
    NA_character_, "picking parasites off a host (cleaner)", NA_character_, "plants/detritus+animals (troph. 2.2-2.79)", 
    "Generalist", "selective plankton feeding", "Planktivorous", "sucking food-containing material", "Benthivorous",
    "variable", "Generalist", NA_character_, NA_character_), by = "FeedingType") |> drop_na(Feeding_mode) |> 
  left_join(rfishbase::load_taxa() |> dplyr::select(SpecCode, Species), by = "SpecCode") |> 
  left_join(tab_Glob |> dplyr::select(Species, Species_FLUXGLOB), by = "Species") |> 
  dplyr::select(- c(SpecCode, FeedingType, Species)) |> relocate(c(Species_FLUXGLOB), .before = Feeding_mode) |> 
  drop_na(Species_FLUXGLOB)

# Merge them all
Guilds <- tab_Glob |> dplyr::select(Species_FLUXGLOB) |> distinct() |>
  left_join(rbind(Feeding_habits, dTraits, MEDFLUX) |> 
              mutate(Feeding_mode = stringr::str_to_sentence(Feeding_mode)) |>
      distinct(Species_FLUXGLOB, Feeding_mode) |> group_by(Species_FLUXGLOB) |>
      summarise(Feeding_mode = paste(sort(unique(Feeding_mode)), collapse = ", "), 
                .groups = "drop"), by = "Species_FLUXGLOB")

# Extract document to fill manually with Experts
writexl::write_xlsx(Guilds, "Outputs/FLUXGLOB/dat_proc/Feeding_mode_to_fill.xlsx")