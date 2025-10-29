#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom") ; library("gstat") ; library("forcats")

## Download data
Env <- read.delim("Data/HAULS_DB_ENV_HINDCAST_COMPLETE_2024.csv", sep = ";")

## Charge from previous scripts
load("Data/HAULS_DB_FPI.RData")
load("Outputs/dat_proc/medits_sf_percentile.Rdata")
load("Outputs/dat_proc/Medit_FunCatch_without_NA.RData")
load("Outputs/dat_proc/Temp_change_slopes.RData")

## Functions
source("Scripts/00_functions_script.R")

Env_FI = merge(Env, HAULS_DB_FPI, by = c("ids", "X", "Y", "YEAR")) |> 
  mutate(GSA = as.character(merged_sf_mult$GSA), MONTH = as.character(MONTH))
merged_sf <- st_join(Slopes, st_as_sf(Env_FI, coords = c("X", "Y"), crs = 4326), join = st_nearest_feature) |> 
  dplyr::filter(Mf_trend_class %in% c("4. significantly positive", "1. significantly negative")) |> 
  mutate(FPI_combined = coalesce(nearest_FPI, FPI_tot)) |> select(-c(nearest_FPI, FPI_tot)) |> 
  mutate(Mf_bin = ifelse(Multifunctionality_slope > 0, 1, 0), abs_botTemp_anom = abs(botTemp_anom))

model = brms::brm(formula = Mf_bin ~ botTemp_anom + DEPTH + chl + FPI_combined + (1|GSA), 
                    data = merged_sf, cores = 3, chains = 3, 
                    family = bernoulli(link = "logit"))
(post_summary <- as.data.frame(fixef(model)))
bayes_R2(model)

ce_FPI <- conditional_effects(model, effects = "FPI_combined")
ce_sbta <- conditional_effects(model, effects = "botTemp_anom")

# Plot
plot(ce_FPI, points = TRUE)
plot(ce_sbta, points = TRUE)
