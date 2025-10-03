#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Hexagonal_grid     <- st_read("Data/Grid_0-1000m_Med.shp")
Guilds             <- readxl::read_xlsx("Data/Guilds_MED.xlsx") |> mutate(SPECIES = str_replace_all(SPECIES, "\u00A0", " "))

## Charge from previous scripts
load("Outputs/dat_proc/Med_all.RData")
load("Outputs/dat_proc/Medit_Temp.RData")
load("Outputs/dat_proc/medits_sf_percentile.Rdata")
load("Outputs/dat_proc/Medit_FunCatch_without_NA.RData")

## Color palette and map
land       <- ne_countries(scale = "medium", returnclass = "sf")
my_palette <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffff",
               "#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")

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
                        filter(n_years >= 10) |> st_drop_geometry() |> ungroup() |> select(HEX_ID) |> 
                        pull(HEX_ID))) |> group_by(HEX_ID) 

## Analyze slopes over year
Slopes <- st_as_sf(medits_with_hexID_over_10years %>%
  pivot_longer(cols = c(community_Gc, community_Fn, community_Fp, Ic_plank, Ic_benthivorous),
    names_to = "metric", values_to = "value") %>% group_by(HEX_ID, metric) %>% do(tidy(lm(value ~ YEAR, data = .))) %>%
  filter(term == "YEAR") %>% select(HEX_ID, metric, slope = estimate) %>%
  pivot_wider(names_from = metric, values_from = slope, names_glue = "{metric}_slope") %>%
  left_join(medits_with_hexID_over_10years %>% st_drop_geometry() %>% group_by(HEX_ID) %>%
      summarise(centroid_lon = mean((left + right) / 2, na.rm = TRUE),
                centroid_lat = mean((top + bottom) / 2, na.rm = TRUE)) %>%
      st_as_sf(coords = c("centroid_lon", "centroid_lat"), crs = st_crs(medits_with_hexID_over_10years))))


ggplot(Slopes) + 
  geom_sf(aes(fill = community_Gc_slope), color = "black", shape = 21) +
  scale_fill_gradientn(colors = my_palette, 
                       name = "Gc slope",
                       na.value = "grey90") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + 
  coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    plot.title      = element_text(size = 20, hjust = 0),
    axis.title      = element_text(size = 18),
    axis.text       = element_text(size = 16),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = "bottom"
  )
hist(Slopes$community_Fp_slope)
