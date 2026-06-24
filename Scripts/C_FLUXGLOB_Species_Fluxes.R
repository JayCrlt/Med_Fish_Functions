#### Setting up          ----
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales") ; library("data.table")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("ggforce") ; library("ggtext")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("purrr") ; library("sf")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("hexbin") ; library("leaflet.extras") ; library("brms")
library("RColorBrewer") ; library("MASS") ; library("stringr")

## Functions
source("Scripts/00_functions_script.R")
source("Scripts/A_FLUXGLOB_Context.R")

## Save data
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB.RData")
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB_Species_list.RData")

# Nutrients database
Nutrients_FishBase_clean = Nutrients_FishBase |> dplyr::select(-Family) |> 
  left_join(load_taxa(), by = "Species") |> group_by(Species, Genus, Subfamily, Family, Order) |> 
  summarise(`C/DM` = mean(`C/DM`), `N/DM` = mean(`N/DM`), `P/DM` = mean(`P/DM`)) |> 
  mutate(Species = gsub(" ", "_", Species))
phy <- fishtree_phylogeny(species = Nutrients_FishBase_clean$Species) 
phy$tip.label <- gsub(" ", "_", phy$tip.label)
nutrients_imputed <- Nutrients_FishBase_clean |>
  impute_trait_phylopars_FG(trait_col = "C/DM", phy = phy) |>
  impute_trait_phylopars_FG(trait_col = "N/DM", phy = phy) |>
  impute_trait_phylopars_FG(trait_col = "P/DM", phy = phy)

#### 7. Compilation  ----
# Needed for the last phylogenetic imputation with body composition
sp. =  unique(c(Species_list$Species, nutrients_imputed$Species)) 
phy_lw              <- fishtree_phylogeny(species = sp.) # 1576 species out of 1946 (+ ~700 single name)
phy_lw$tip.label    <- gsub(" ", "_", phy_lw$tip.label)

###### 7.1 Clean Species info  ----
tab_Glob <- Species_list |> mutate(SpecCode = as.integer(SpecCode)) |> 
  left_join(Nutrients_FishBase |> dplyr::select(-Family) |> left_join(load_taxa(), by = "Species") |> 
              group_by(Species, Genus, Subfamily, Family, Order) |> 
              summarise(`C/DM` = mean(`C/DM`), `N/DM` = mean(`N/DM`), `P/DM` = mean(`P/DM`)), 
            by = c("Species", "Genus", "Subfamily", "Family", "Order")) |> 
  left_join(rfishbase::popqb(), by = "SpecCode") |> 
  dplyr::select(Species, Genus, Subfamily, Family, Order, `C/DM`, `N/DM`, `P/DM`, PopQB, FoodType) |> 
  rename(Qc = `C/DM`, Qn = `N/DM`, Qp = `P/DM`) |> 
  
  ###### 7.2 Add CNP body composition ----
bind_rows(nutrients_imputed |> dplyr::select(Species, Genus, Subfamily, Family, Order, CDM, NDM, PDM) |> 
            rename(Qc = CDM, Qn = NDM, Qp = PDM) |> mutate(Dataset = "Schiettekatte_2021")) |> 
  filter(!is.na(Species)) |> 
  mutate(Species = gsub(" ", "_", Species)) |> 
  impute_trait_phylopars_FG(trait_col = "Qc", phy = phy_lw) |>
  impute_trait_phylopars_FG(trait_col = "Qn", phy = phy_lw) |>
  impute_trait_phylopars_FG(trait_col = "Qp", phy = phy_lw) |> 
  dplyr::select(-c(Qc_type, Qn_type, Qp_type)) |> relocate(Dataset, .after = Qp) |> group_by(Species) |>
  slice_max(order_by = (Dataset == "FLUXGLOB"), with_ties = FALSE) |> ungroup() |> 
  fill_trait_hierarchy_FG("Qc") |> 
  fill_trait_hierarchy_FG("Qn") |> 
  fill_trait_hierarchy_FG("Qp") |> 
  
  ###### 7.3 Add CNP diet             ----
mutate(Species = gsub("_", " ", Species)) |> 
  left_join(rfishbase::load_taxa(), by = c("Species", "Genus", "Subfamily", "Family", "Order")) |> 
  left_join(rfishbase::ecology(), by = "SpecCode", relationship = "many-to-many") |> 
  dplyr::select(c(1:9,51)) |> mutate(Species = gsub(" ", "_", Species)) |> 
  impute_trait_phylopars_FG(trait_col = "FoodTroph", phy = phy_lw) |> 
  fill_trait_hierarchy_FG("FoodTroph") |> dplyr::select(-FoodTroph_type) |> 
  mutate(Dc = predict(model_C, newdata = pick(everything())),
         Dn = predict(model_N, newdata = pick(everything())),
         Dp = predict(model_P, newdata = pick(everything()))) |> 
  relocate(c(FoodTroph, Dataset), .after = Dp) |> mutate(Qc = Qc * 100, Qn = Qn * 100, Qp = Qp * 100) |> 
  dplyr::filter(Dataset == "B_Useful_MED") |> rename(h = FoodTroph) |> 
  
  ###### 7.4 Add constantes           ----
mutate(ac = 0.8, an = 0.8, ap = 0.7, F0Nz = 3.7e-03, F0Pz = 3.7e-04) |> 
  relocate(Dataset, .after = F0Pz) |> 
  left_join(mdw_FishBase_imputed |> mutate(Species = gsub(" ", "_", Species))) |> 
  fill_trait_hierarchy("mdw") |> 
  relocate(Dataset, .after = mdw) |> 
  
  ###### 7.5 Add growth parameters    ----
mutate(lwa = NA, lwb = NA, linf = NA, K = NA, t = NA) |> 
  full_join(lw_growth |> dplyr::select(Species, lwa, lwb,linf, K, t) |> 
              mutate(Dataset = "FishBase") |> rename(t0 = t) |> 
              relocate(Dataset, .after = t0)) |> 
  relocate(Dataset, .after = t0) |> group_by(Species) |> 
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop")  |> 
  fill_trait_hierarchy("lwa") |> fill_trait_hierarchy("lwb") |> fill_trait_hierarchy("linf") |> 
  fill_trait_hierarchy("K") |> fill_trait_hierarchy("t0") |> 
  
  ###### 7.6 Add Caudal fin          ----
mutate(r = NA) |> 
  full_join(Caudal_fin_FB |> dplyr::select(Species, AspectRatio) |> 
              rename(r = AspectRatio) |> 
              mutate(Dataset = "FishBase") |> 
              relocate(Dataset, .after = r)) |> 
  relocate(Dataset, .after = r) |> group_by(Species) |> 
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop")  |> 
  fill_trait_hierarchy("r") |> 
  
  ###### 7.6 Add Metabolism data     ----
# Workflow: get SMR and MMR thanks to relationships obtained in section 6.0. 
# Then, use: f0 = SMR_gC_day / Weight^alpha * exp(E * (1 / (Temp_ref + 273.15) - 1 / (Temp_ref + 273.15)) / 8.617e-5))
# And theta  = (SMR_gC_day + MMR_gC_day) / (2 * SMR_gC_day)
mutate(alpha  = round(fixef(SMR_Weight_relationship)["log_weight", "Estimate"], 3)) |> 
  dplyr::filter(Dataset == "B_Useful_MED") |> 
  relocate(Dataset, .after = alpha) |> dplyr::select(-t) |> 
  
  ##### 7.7 Format data to match obs ----
mutate(Species = case_when(
  Species == "Chelidonichthys lastoviza" ~ "Trigloporus lastoviza",
  Species == "Aetomylaeus bovinus"       ~ "Pteromylaeus bovinus",
  Species == "Chelon ramada"             ~ "Liza ramada",
  Species == "Chelon auratus"            ~ "Liza aurata",
  Species == "Chelon saliens"            ~ "Liza saliens",
  TRUE ~ Species))

#### 8. Export the data  ----
## Data
save(All_Med, file = "Outputs/dat_proc/All_Med.RData")