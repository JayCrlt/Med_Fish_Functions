#### Setting up          ----
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales") ; library("data.table")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("ggforce") ; library("ggtext")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("purrr") ; library("sf")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("hexbin") ; library("leaflet.extras") ; library("brms")
library("RColorBrewer") ; library("MASS") ; library("stringr") ; library("terra") ; library("reticulate")

## Functions
source("Scripts/00_functions_script.R")
source("Scripts/A_FLUXGLOB_Context.R")
source("Scripts/B_FLUXGLOB_Cleaning_Species.R")

## Load data
load("Outputs/dat_proc/model_C.RData")
load("Outputs/dat_proc/model_N.RData")
load("Outputs/dat_proc/model_P.RData")
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB.RData")
load("Outputs/dat_proc/SMR_Weight_relationship.RData")
load("Outputs/dat_proc/MMR_Weight_relationship.RData")
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB_Species_list.RData")

## Nutrients, growth database and Caudal fin database
# Nutrients
Nutrients_FishBase_clean = Nutrients_FishBase |> dplyr::select(-Family) |> 
  left_join(load_taxa(), by = "Species") |> group_by(Species, Genus, Subfamily, Family, Order) |> 
  summarise(CDM = mean(`C/DM`), NDM = mean(`N/DM`), PDM = mean(`P/DM`)) |> 
  mutate(Species = gsub(" ", "_", Species))
phy <- fishtree_phylogeny(species = Nutrients_FishBase_clean$Species) 
phy$tip.label <- gsub(" ", "_", phy$tip.label)
nutrients_imputed <- Nutrients_FishBase_clean |>
  impute_trait_phylopars_FG(trait_col = "CDM", phy = phy) |>
  impute_trait_phylopars_FG(trait_col = "NDM", phy = phy) |>
  impute_trait_phylopars_FG(trait_col = "PDM", phy = phy)

# Growth
lw_growth <- rfishbase::length_weight() |> group_by(SpecCode) |>
  summarise(lwa = mean(a, na.rm = T), lwb = mean(b, na.rm = T)) |>
  left_join(rfishbase::popgrowth() |> group_by(SpecCode) |>
              summarise(linf = mean(Loo, na.rm = T), K = mean(K, na.rm = T), t0 = mean(to, na.rm = T)), by = "SpecCode") |>
  left_join(rfishbase::load_taxa() |> dplyr::select(SpecCode, Species, Genus, Subfamily, Family, Order),  by = "SpecCode") |> 
  dplyr::select(Species, everything()) |> mutate(t0 = ifelse(is.nan(t0), NA, t0)) |> 
  dplyr::select(-SpecCode)

lw_growth$Species <- gsub(" ", "_", lw_growth$Species)
phy_lw            <- fishtree_phylogeny(species = lw_growth$Species)
phy_lw$tip.label  <- gsub(" ", "_", phy_lw$tip.label)
lw_growth         <- lw_growth |> filter(!is.na(Species))

lw_growth         <- lw_growth |>
  impute_trait_phylopars_FG(trait_col = "lwa", phy = phy_lw) |>
  impute_trait_phylopars_FG(trait_col = "lwb", phy = phy_lw) |>
  impute_trait_phylopars_FG(trait_col = "linf", phy = phy_lw) |> 
  impute_trait_phylopars_FG(trait_col = "K", phy = phy_lw) |>
  impute_trait_phylopars_FG(trait_col = "t0", phy = phy_lw) 

# Caudal fin
Caudal_fin_FB <- rfishbase::morphometrics() |> data.frame() |> 
  dplyr::select(SpecCode, AspectRatio) |> 
  left_join(rfishbase::load_taxa(), by = "SpecCode") |> 
  dplyr::select(Species, Genus, Subfamily, Family, Order, AspectRatio) |> 
  group_by(Species, Genus, Subfamily, Family, Order) |> summarise(AspectRatio = mean(AspectRatio, na.rm = TRUE)) |>
  mutate(AspectRatio = ifelse(is.nan(AspectRatio), NA, AspectRatio)) |> ungroup()

Caudal_fin_FB$Species <- gsub(" ", "_", Caudal_fin_FB$Species)
phy_caudal <- fishtree_phylogeny(species = Caudal_fin_FB$Species)
phy_caudal$tip.label <- gsub(" ", "_", phy_caudal$tip.label)
caudal_imputed <- Caudal_fin_FB |> filter(!is.na(Species)) |> 
  impute_trait_phylopars_FG(trait_col = "AspectRatio", phy = phy_caudal)

#### Compilation  ----
# Needed for the last phylogenetic imputation with body composition
sp. =  unique(c(Species_list$Species, nutrients_imputed$Species)) 
phy_lw              <- fishtree_phylogeny(species = sp.) # 1576 species out of 1946 (+ ~700 single name)
phy_lw$tip.label    <- gsub(" ", "_", phy_lw$tip.label)

###### 1. Clean Species info  ----
tab_Glob <- Species_list |> mutate(SpecCode = as.integer(SpecCode)) |> 
  left_join(Nutrients_FishBase |> dplyr::select(-Family) |> left_join(load_taxa(), by = "Species") |> 
              group_by(Species, Genus, Subfamily, Family, Order) |> 
              summarise(`C/DM` = mean(`C/DM`), `N/DM` = mean(`N/DM`), `P/DM` = mean(`P/DM`)), 
            by = c("Species", "Genus", "Subfamily", "Family", "Order")) |> 
  left_join(rfishbase::popqb(), by = "SpecCode") |> 
  dplyr::select(Species, Genus, Subfamily, Family, Order, `C/DM`, `N/DM`, `P/DM`, PopQB, FoodType) |> 
  rename(Qc = `C/DM`, Qn = `N/DM`, Qp = `P/DM`) |> 
  
  ###### 2. Add CNP body composition ----
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
  
  ###### 3. Add CNP diet             ----
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
  rename(h = FoodTroph) |> 
  
  ###### 4. Add constantes           ----
  mutate(ac = 0.8, an = 0.8, ap = 0.7, F0Nz = 3.7e-03, F0Pz = 3.7e-04) |> 
  relocate(Dataset, .after = F0Pz) |> 
  left_join(mdw_FishBase_imputed |> mutate(Species = gsub(" ", "_", Species))) |> 
  fill_trait_hierarchy_FG("mdw") |> 
  relocate(Dataset, .after = mdw) |> 
  
  ###### 5. Add growth parameters    ----
  mutate(lwa = NA, lwb = NA, linf = NA, K = NA, t0 = NA) |> 
  full_join(lw_growth |> dplyr::select(Species, lwa, lwb,linf, K, t0) |> 
              mutate(Dataset = "FishBase") |> 
              relocate(Dataset, .after = t0)) |> 
  relocate(Dataset, .after = t0) |> group_by(Species) |> 
  fill_trait_hierarchy_FG("lwa") |> fill_trait_hierarchy_FG("lwb") |> fill_trait_hierarchy_FG("linf") |> 
  fill_trait_hierarchy_FG("K") |> fill_trait_hierarchy_FG("t0") |> ungroup() |> 
  
  ###### 6. Add Caudal fin           ----
  mutate(r = NA) |> 
  full_join(Caudal_fin_FB |> dplyr::select(Species, AspectRatio) |> 
              rename(r = AspectRatio) |> 
              mutate(Dataset = "FishBase") |> 
              relocate(Dataset, .after = r)) |> 
  relocate(Dataset, .after = r) |> group_by(Species) |> 
  fill_trait_hierarchy_FG("r") |> ungroup() |> 
  
  ###### 7.  Add Metabolism data     ----
# Workflow: get SMR and MMR thanks to relationships obtained in section 6.0. 
# Then, use: f0 = SMR_gC_day / Weight^alpha * exp(E * (1 / (Temp_ref + 273.15) - 1 / (Temp_ref + 273.15)) / 8.617e-5))
# And theta  = (SMR_gC_day + MMR_gC_day) / (2 * SMR_gC_day)
  mutate(alpha  = round(fixef(SMR_Weight_relationship)["log_weight", "Estimate"], 3)) |> 
  relocate(Dataset, .after = alpha) |> 
  
  ##### 8.  Format data to match obs ----
  mutate(Species = gsub("_", " ", Species)) |> dplyr::filter(Species %in% FLUXGLOB_cor$Species) |> 
  dplyr::select(-c(Subfamily, Dataset)) |> ungroup() |> group_by(Species) |> 
  summarise(across(everything(), ~ coalesce(.x[1], .x[which(!is.na(.x))[1]])), .groups = "drop") |> 
  dplyr::select(-c("Genus", "Family", "Order")) |> left_join(Species_list |> dplyr::select(-c(SpecCode))) |> 
  relocate(c(Genus, Subfamily, Family, Order), .after = Species) |> relocate(c(Species_FLUXGLOB), .before = Species) |> 
  dplyr::select(-c("Class", "SuperClass")) |> drop_na(Species) # Remove Selachi and Perciformes without info from df

#### Export the data  ----
save(tab_Glob, file = "Outputs/FLUXGLOB/dat_proc/All_Species_FLUXGLOB.RData")