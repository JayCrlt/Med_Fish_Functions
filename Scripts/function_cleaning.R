#### Setting up          ----
rm(list = ls())
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Medits_total       <- read.delim("Data/TATB_WMED_1999-2021_clean.csv", sep = ";") |> 
  left_join(Temp_WMED <- read.delim("Data/HAULS_DB_ENV_HINDCAST_COMPLETE_2024.csv", sep = ";") |> 
    select(-id) |> rename(id = ids) |> select(c(id,botTemp_anom)) |> rename(Temp_anomaly = botTemp_anom))
Hexagonal_grid     <- st_read("Data/Grid 0-1000m_WMED/Grid 0-1000m_WMED.shp")

## Charge from previous scripts
load("Outputs/dat_proc/Western_Med.RData")

## Functions
source("Scripts/functions_script.R")

# add spatial hex info
medits_sf               <- st_as_sf(data.frame(Longitude = as.numeric(gsub(",", ".", Medits_total$MEAN_LONGITUDE_DEC)), 
                                               Latitude = as.numeric(gsub(",", ".", Medits_total$MEAN_LATITUDE_DEC))), 
                                    coords = c("Longitude", "Latitude"), crs = 4326)
medits_sf_hex           <- st_join(medits_sf, Hexagonal_grid["grid.id"], left = TRUE)
Medits_total$HEX_ID     <- medits_sf_hex$grid.id
Medits_total$MASS_IND   <- Medits_total$TOTAL_WEIGHT_IN_THE_HAUL / Medits_total$TOTAL_NUMBER_IN_THE_HAUL

# Clean data according to biomass vs species
Medits_total            <- Medits_total |> dplyr::filter(TOTAL_WEIGHT_IN_THE_HAUL != 0 & TOTAL_NUMBER_IN_THE_HAUL != 0) |> 
  left_join(sp_code_list |> rename(MEDITS_CODE = spp.code)) |>
  filter(catfau %in% c("Aa", "Ae", "Ao")) 

# Some species are unexpectedly heavy...
Western_Med$Species  <- gsub(" ", "_", Western_Med$Species)
phy_lw               <- fishtree_phylogeny(species = unique(Western_Med$Species))
phy_lw$tip.label     <- gsub(" ", "_", phy_lw$tip.label)

Western_Med = Western_Med |> 
  mutate(Species = gsub("_", " ", Species)) |> 
  left_join((rfishbase::load_taxa() |> select(Species, SpecCode)), by = "Species") |> 
  mutate(Genus = ifelse(grepl("spp\\.", Species), word(Species, 1), Genus)) |> 
  left_join((rfishbase::popchar() |> select(SpecCode, Wmax)), by = "SpecCode", relationship = "many-to-many") |> 
  group_by(Species) |> slice_max(order_by = Wmax, n = 1, with_ties = FALSE) |> ungroup() |> 
  left_join((rfishbase::popgrowth() |> select(SpecCode, Winfinity)), by = "SpecCode", relationship = "many-to-many") |> 
  group_by(Species) |> slice_max(order_by = Winfinity, n = 1, with_ties = FALSE) |> ungroup() |> 
  impute_trait_phylopars(trait_col = "Wmax", phy = phy_lw) |> fill_trait_hierarchy("Wmax") |> 
  impute_trait_phylopars(trait_col = "Winfinity", phy = phy_lw) |> fill_trait_hierarchy("Winfinity") |> 
  select(-c(SpecCode, Wmax_type, Winfinity_type)) |> relocate(Dataset, .after = Winfinity) 

Medits_total            <- Medits_total |> 
  mutate(sci.name = case_when(
    sci.name == "Trigloporus lastoviza" ~ "Chelidonichthys lastoviza",
    sci.name == "Pteromylaeus bovinus"  ~ "Aetomylaeus bovinus",
    sci.name == "Liza ramada"           ~ "Chelon ramada",
    sci.name == "Liza aurata"           ~ "Chelon auratus",
    sci.name == "Liza saliens"          ~ "Chelon saliens",
    TRUE ~ sci.name)) |> left_join(Western_Med |> rename(sci.name = Species)) |>
  filter(!grepl("^NO\\s+", sci.name))

# Clean errors according to maximum weight
Medit_Western_FunCatch <- Medits_total |> 
  mutate(NegativeDifference = (if_else(Winfinity > Wmax, Winfinity, Wmax) - MASS_IND) < -100) |> # 100g treshold
  dplyr::filter(NegativeDifference == FALSE) |>  # -72 rows 
  dplyr::select(c(3:4, 9:10, 12, 17, 22, 24, 27:29, 32:33, 37:41, 47:50, 52:72)) |> 
  relocate(c(HEX_ID, MEDITS_CODE), .before = COUNTRY) |> 
  rename(SPECIES = sci.name, GENUS = Genus, FAMILY = Family) |> 
  mutate(SMR = 0.002 * MASS_IND^alpha, MMR = 0.006 * MASS_IND^0.779,    
         TEMP = rowMeans(across(c(BOTTOM_TEMPERATURE_BEGINNING, BOTTOM_TEMPERATURE_END), ~ {
             temp_num <- as.numeric(gsub(",", ".", .)); ifelse(temp_num > 100, temp_num / 10, temp_num)}), na.rm = TRUE),
         f0 = SMR / MASS_IND^alpha * exp(-0.138 * (1 / (TEMP + 273.15) - 1 / (TEMP + Temp_anomaly + 273.15)) / 8.617e-5),
         f0 = if_else(is.nan(f0), mean(f0[!is.nan(f0)], na.rm = TRUE), f0),
         theta  = (SMR + MMR) / (2 * SMR)) |> 
  select(-c(SMR, MMR, BOTTOM_TEMPERATURE_BEGINNING, BOTTOM_TEMPERATURE_END, Temp_anomaly)) |> 
  relocate(TEMP, .after = DISTANCE) |> mutate(mdw = 0.3) |> relocate(mdw, .after = lwb)

##### NEED TO WORK ON mdw PARAMETER
# https://www.fishbase.se/summary/biomass.php?id=69 

# To process the data, split into hex cell and years
Medit_Western_FunCatch_split <- Medit_Western_FunCatch |> group_by(HEX_ID, YEAR) |> group_split() # 11172 docs