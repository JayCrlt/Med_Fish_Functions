#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Hexagonal_grid     <- st_read("Data/Grid_0-1000m_Med.shp")
Guilds             <- readxl::read_xlsx("Data/Guilds_MED.xlsx") |> mutate(SPECIES = str_replace_all(SPECIES, "\u00A0", " "))

## Charge from previous scripts
load("Outputs/dat_proc/Med_all.RData")
load("Outputs/dat_proc/Medit_Temp.RData")
load("Outputs/dat_proc/medits_sf_percentile.Rdata")
load("Outputs/dat_proc/Medit_FunCatch_without_NA.RData")

## Color palette
c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffff","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")

## Functions
source("Scripts/00_functions_script.R")

## Check visually hexagon grid with data
Hexagonal_grid_with_points <- st_join(st_transform(Hexagonal_grid, st_crs(medits_sf_percentile)), 
                                      medits_sf_percentile, join = st_intersects) |> drop_na(YEAR) |> 
  ggplot() + geom_sf(fill = "red", alpha = 0.4, color = "black") + theme_minimal() +
  geom_sf(data = ne_countries(scale = "medium", returnclass = "sf"), fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

## Attribute an ID from Hexagonal_grid to dataset
medits_with_hexID <- st_join(medits_sf_percentile, st_transform(Hexagonal_grid, st_crs(medits_sf_percentile)), 
  join = st_intersects) |> drop_na(id) |> rename(HEX_ID = id) |> relocate(HEX_ID, .after = YEAR)

## Keep only HEX cells with at least 10 years of data
medits_with_hexID_over_10years <- medits_with_hexID |> 
  filter(HEX_ID %in% (medits_with_hexID |> group_by(HEX_ID) |> summarise(n_years = n_distinct(YEAR)) |> 
                        filter(n_years >= 10) |> st_drop_geometry() |> ungroup() |> select(HEX_ID) |> pull(HEX_ID))) |> 
  group_by(HEX_ID) |> 
  # Define outliers for Gc
  mutate(Q1 = quantile(community_Gc, 0.25, na.rm = T), Q3 = quantile(community_Gc, 0.75, na.rm = T), IQR = Q3 - Q1,
  outlier_Gc = ifelse(community_Gc < (Q1 - 1.5 * IQR) | community_Gc > (Q3 + 1.5 * IQR), "outlier", NA)) |> 
  select(-Q1, -Q3, -IQR) |> ungroup() |> 
  # Define outliers for Fn
  mutate(Q1 = quantile(community_Fn, 0.25, na.rm = T), Q3 = quantile(community_Fn, 0.75, na.rm = T), IQR = Q3 - Q1,
  outlier_Fn = ifelse(community_Fn < (Q1 - 1.5 * IQR) | community_Fn > (Q3 + 1.5 * IQR), "outlier", NA)) |> 
  select(-Q1, -Q3, -IQR) |> ungroup() |> 
  # Define outliers for Fp
  mutate(Q1 = quantile(community_Fp, 0.25, na.rm = T), Q3 = quantile(community_Fp, 0.75, na.rm = T), IQR = Q3 - Q1,
  outlier_Fp = ifelse(community_Fp < (Q1 - 1.5 * IQR) | community_Fp > (Q3 + 1.5 * IQR), "outlier", NA)) |> 
  select(-Q1, -Q3, -IQR) |> ungroup() |>  
  # Define outliers for Plank
  mutate(Q1 = quantile(Ic_plank, 0.25, na.rm = T), Q3 = quantile(Ic_plank, 0.75, na.rm = T), IQR = Q3 - Q1,
  outlier_plank = ifelse(Ic_plank < (Q1 - 1.5 * IQR) | Ic_plank > (Q3 + 1.5 * IQR), "outlier", NA)) |> 
  select(-Q1, -Q3, -IQR) |> ungroup() |> 
  # Define outliers for Benth
  mutate(Q1 = quantile(Ic_benthivorous, 0.25, na.rm = T), Q3 = quantile(Ic_benthivorous, 0.75, na.rm = T), IQR = Q3 - Q1,
  outlier_benthivorous = ifelse(Ic_benth < (Q1 - 1.5 * IQR) | Ic_benth > (Q3 + 1.5 * IQR), "outlier", NA)) |> 
  select(-Q1, -Q3, -IQR) |> ungroup()

hist(log(medits_sf_percentile$Ic_benthivorous))

## Analyze Gc slopes over year
ggplot(data = medits_with_hexID_over_10years |> dplyr::filter(HEX_ID %in% c(1057, 1177)), 
       aes(y = community_Gc, x = YEAR, group = HEX_ID)) + geom_line() +
  geom_point(aes(fill = outlier_Gc), shape = 21, size = 3)
