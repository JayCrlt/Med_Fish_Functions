rm(list = ls())
library("rfishbase") ; library("phylosem") ; library("tidyverse")

## What is a function? A function is a rate. (Jax)
## Then, the idea would be to define three core functions: i) Biomass production through Growth rates;
## ii) Predation rate through  Q/B * Biomass; iii) nutrient cycling through C,N,P/DM + FishFlux
## All those data are available in FishBase, but maybe not directly from RFishBase, then, we would have to compile 
## The data by ourselve

# Raw data
sp_code_list  <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Medits_total  <- read.delim("Data/TATB_WMED_1999-2021_clean.csv", sep = ";")

# Merge Species_code with Scientific name
Medits_total = Medits_total |> mutate(spp.code = paste(GENUS, SPECIES, sep = "")) |> 
  left_join(sp_code_list, by = "spp.code") |> 
  filter(!str_detect(sci.name, "^NO\\s*")) |> 
  filter(catfau %in% c("Aa", "Ae", "Ao")) # Work only with fish first

# How many species in the Med Sea to play with?
length(unique(Medits_total$sci.name)) # 361 fish out of 1246 species in total

# Some Viz
(Biomass_GSA_year <- Medits_total |> group_by(GSA, YEAR) |> 
  summarise(sum(TOTAL_WEIGHT_IN_THE_HAUL)) |> 
  ggplot(aes(y = `sum(TOTAL_WEIGHT_IN_THE_HAUL)`, x = YEAR)) +
  geom_line() + facet_wrap(~GSA))
