#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("rphylopic")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Hexagonal_grid     <- st_read("Data/Grid_0-1000m_Med.shp")
Guilds             <- readxl::read_xlsx("Data/Guilds_MED.xlsx") |> mutate(SPECIES = str_replace_all(SPECIES, "\u00A0", " "))

## Charge from previous scripts
load("Outputs/dat_proc/Med_all.RData")
load("Outputs/dat_proc/Medit_Temp.RData")

## Functions
source("Scripts/00_functions_script.R")

# add spatial hex info
medits_sf               <- st_as_sf(data.frame(Longitude = as.numeric(gsub(",", ".", MEDIT_TEMP$MEAN_LONGITUDE_DEC)), 
                                               Latitude  = as.numeric(gsub(",", ".", MEDIT_TEMP$MEAN_LATITUDE_DEC)),
                                               ID        = seq_len(nrow(MEDIT_TEMP))), 
                                    coords = c("Longitude", "Latitude"), crs = 4326)
medits_sf_hex           <- st_join(medits_sf, Hexagonal_grid["grid_id"], left = TRUE)
medits_sf_hex_unique    <- medits_sf_hex |> st_drop_geometry() |> group_by(ID) |> summarise(grid_id = first(grid_id)) |> ungroup()
MEDIT_TEMP$HEX_ID       <- medits_sf_hex_unique$grid_id
MEDIT_TEMP$MASS_IND     <- MEDIT_TEMP$TOTAL_WEIGHT_IN_THE_HAUL / MEDIT_TEMP$TOTAL_NUMBER_IN_THE_HAUL

# Clean data according to biomass vs species
Medits_total            <- MEDIT_TEMP |> dplyr::filter(TOTAL_WEIGHT_IN_THE_HAUL != 0 & TOTAL_NUMBER_IN_THE_HAUL != 0) |> 
  left_join(sp_code_list |> rename(MEDITS_CODE = spp.code)) |>
  filter(catfau %in% c("Aa", "Ae", "Ao")) 

# Some species are unexpectedly heavy...
Med_all$Species  <- gsub(" ", "_", Med_all$Species)
phy_lw               <- fishtree_phylogeny(species = unique(Med_all$Species))
phy_lw$tip.label     <- gsub(" ", "_", phy_lw$tip.label)

Med_all = Med_all |> 
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
    TRUE ~ sci.name)) |> left_join(Med_all |> rename(sci.name = Species)) |>
  filter(!grepl("^NO\\s+", sci.name))

# Clean errors according to maximum weight
Medit_FunCatch <- Medits_total |> 
  mutate(NegativeDifference = (if_else(Winfinity > Wmax, Winfinity, Wmax) - MASS_IND) < -100) |> # 100g treshold
  dplyr::filter(NegativeDifference == FALSE) |>  # -72 rows 
  dplyr::select(c(3:4, 9:12, 17, 22, 24:29, 32:33, 37:41, 47:51, 53:73)) |> 
  relocate(c(HEX_ID, MEDITS_CODE), .before = COUNTRY) |> 
  rename(SPECIES = sci.name, GENUS = Genus, FAMILY = Family) |> 
  mutate(SMR = 0.002 * MASS_IND^0.836, MMR = 0.006 * MASS_IND^0.779,    
         f0 = SMR / MASS_IND^0.836 * exp(-0.138 * (1 / (TEMP + 273.15) - 1 / (TEMP + Temp_anomaly + 273.15)) / 8.617e-5),
         f0 = if_else(is.nan(f0), mean(f0[!is.nan(f0)], na.rm = TRUE), f0),
         theta  = (SMR + MMR) / (2 * SMR),
         HAUL_AREA = as.numeric(DISTANCE) * WING_OPENING,
         HAUL_VOLUME = as.numeric(DISTANCE) * WING_OPENING * VERTICAL_OPENING) |> 
  select(-c(SMR, MMR, BOTTOM_TEMPERATURE_BEGINNING, BOTTOM_TEMPERATURE_END, Temp_anomaly, WING_OPENING, VERTICAL_OPENING)) |> 
  relocate(c(HAUL_AREA, HAUL_VOLUME, TEMP), .after = DISTANCE) |> relocate(mdw, .after = lwb) |> 
  mutate(TL = (MASS_IND / lwa)^(1 / lwb))

## FishFlux loop workflow
# Prepare storing vectors
model         <- list()
param_dataset <- list()
# Prepare the dataset
Medit_FunCatch = Medit_FunCatch |> 
  mutate(Fn_mean = NA, Fn_Q1 = NA, Fn_Q3 = NA, Fp_mean = NA, Fp_Q1 = NA, Fp_Q3 = NA, 
         Gc_mean = NA, Gc_Q1 = NA, Gc_Q3 = NA, Ic_mean = NA, Ic_Q1 = NA, Ic_Q3 = NA)

# Initiate the FishFlux formatting
for(i in 1:length(Medit_FunCatch$HEX_ID)){
  param_dataset[[i]] <- list(
    ac_m    = Medit_FunCatch$ac[i],
    ac_sd   = Medit_FunCatch$ac[i] * 0.05,
    an_m    = Medit_FunCatch$an[i],
    an_sd   = Medit_FunCatch$an[i] * 0.05,
    ap_m    = Medit_FunCatch$ap[i],
    ap_sd   = Medit_FunCatch$ap[i] * 0.05,
    r_m     = Medit_FunCatch$r[i],
    f0_m    = Medit_FunCatch$f0[i],
    theta_m = Medit_FunCatch$theta[i],
    Dc_m    = Medit_FunCatch$Dc[i],
    Dc_sd   = Medit_FunCatch$Dc[i] * 0.05,
    Dn_m    = Medit_FunCatch$Dn[i],
    Dn_sd   = Medit_FunCatch$Dn[i] * 0.05,
    Dp_m    = Medit_FunCatch$Dp[i],
    Dp_sd   = Medit_FunCatch$Dp[i] * 0.05,
    k_m     = Medit_FunCatch$K[i],
    linf_m  = Medit_FunCatch$linf[i],
    lwa_m   = Medit_FunCatch$lwa[i],
    lwb_m   = Medit_FunCatch$lwb[i],
    Qc_m    = Medit_FunCatch$Qc[i],
    Qc_sd   = Medit_FunCatch$Qc[i] * 0.05,
    Qn_m    = Medit_FunCatch$Qn[i],
    Qn_sd   = Medit_FunCatch$Qn[i] * 0.05,
    Qp_m    = Medit_FunCatch$Qp[i],
    Qp_sd   = Medit_FunCatch$Qp[i] * 0.05,
    t0_m    = Medit_FunCatch$t0[i],
    h_m     = Medit_FunCatch$h[i],
    F0nz_m  = Medit_FunCatch$F0Nz[i],
    F0pz_m  = Medit_FunCatch$F0Pz[i],
    mdw_m   = Medit_FunCatch$mdw[i],
    v_m     = Medit_FunCatch$TEMP[i],
    alpha_m = 0.836)}

# Compile data from FishFlux
pb <- txtProgressBar(min = 0, max = length(Medit_FunCatch$HEX_ID), style = 3)
for(i in 1: length(Medit_FunCatch$HEX_ID)){
  tryCatch({
    invisible(capture.output({
  model[[i]] <- extract(cnp_model_mcmc(TL = Medit_FunCatch$TL[i], param = param_dataset[[i]]), c("Fn","Fp", "Gc", "Ic"))}))
  Medit_FunCatch$Fn_mean[i] <- model[[i]]$Fn_mean
  Medit_FunCatch$Fn_Q1[i]   <- model[[i]]$`Fn_2.5%`
  Medit_FunCatch$Fn_Q3[i]   <- model[[i]]$`Fn_97.5%`
  Medit_FunCatch$Fp_mean[i] <- model[[i]]$Fp_mean
  Medit_FunCatch$Fp_Q1[i]   <- model[[i]]$`Fp_2.5%`
  Medit_FunCatch$Fp_Q3[i]   <- model[[i]]$`Fp_97.5%`
  Medit_FunCatch$Gc_mean[i] <- model[[i]]$Gc_mean
  Medit_FunCatch$Gc_Q1[i]   <- model[[i]]$`Gc_2.5%`
  Medit_FunCatch$Gc_Q3[i]   <- model[[i]]$`Gc_97.5%`
  Medit_FunCatch$Ic_mean[i] <- model[[i]]$Ic_mean
  Medit_FunCatch$Ic_Q1[i]   <- model[[i]]$`Ic_2.5%`
  Medit_FunCatch$Ic_Q3[i]   <- model[[i]]$`Ic_97.5%`
  }, error = function(e){
    message("Error in row ", i)
    Medit_FunCatch$Fn_mean[i] <- NA
    Medit_FunCatch$Fn_Q1[i]   <- NA
    Medit_FunCatch$Fn_Q3[i]   <- NA
    Medit_FunCatch$Fp_mean[i] <- NA
    Medit_FunCatch$Fp_Q1[i]   <- NA
    Medit_FunCatch$Fp_Q3[i]   <- NA
    Medit_FunCatch$Gc_mean[i] <- NA
    Medit_FunCatch$Gc_Q1[i]   <- NA
    Medit_FunCatch$Gc_Q3[i]   <- NA
    Medit_FunCatch$Ic_mean[i] <- NA
    Medit_FunCatch$Ic_Q1[i]   <- NA
    Medit_FunCatch$Ic_Q3[i]   <- NA
  })
  setTxtProgressBar(pb, i)}

# There was significant numbers of errors due to missing observed temperature, it has been fixed thanks to R_script_02
Medit_FunCatch_without_NA = Medit_FunCatch |> 
  anti_join(Medit_FunCatch |> group_by(HEX_ID, HAUL_NUMBER, YEAR, MONTH, HAULING_TIME) |> 
              summarise(NaN_count = sum(is.nan(TEMP)), Total_count = n(), .groups = "drop") |> 
              filter(NaN_count > 0) |> select(HEX_ID, YEAR), by = c("HEX_ID", "YEAR")) # Remove 501533 - 501403 = -130 obs (< -0.001%) 

# Here added guilds
Medit_FunCatch_without_NA <- Medit_FunCatch_without_NA |> full_join(Guilds) |> 
  drop_na(Fn_mean, Fp_mean, Gc_mean, Ic_mean)

# Community computation at the individual level
Medit_FunCatch_without_NA_community = Medit_FunCatch_without_NA |> 
  mutate(SWEPT_AREA = as.numeric(gsub(",", ".", SWEPT_AREA)),
         Biomass = TOTAL_WEIGHT_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 / 24),
         community_Fn = Fn_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 / 24),
         community_Fp = Fp_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 / 24),
         community_Gc = Gc_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 / 24)) |> 
  group_by(HEX_ID, HAUL_NUMBER, YEAR, MONTH, HAULING_TIME, MEAN_LONGITUDE_DEC, MEAN_LATITUDE_DEC) |> 
  summarise(Biomass = sum(Biomass), community_Fn = sum(community_Fn), community_Fp = sum(community_Fp),
            community_Gc = sum(community_Gc), .groups = "drop") |> 
  mutate(top1_Biomass      = rank(-Biomass, ties.method = "first") <= 127, #12,666 rows
         top1_community_Fn = rank(-community_Fn, ties.method = "first") <= 127,
         top1_community_Fp = rank(-community_Fp, ties.method = "first") <= 127,
         top1_community_Gc = rank(-community_Gc, ties.method = "first") <= 127)

# Community computation at the trophic guild level
Medit_FunCatch_without_NA_trophic = Medit_FunCatch_without_NA |> 
  mutate(SWEPT_AREA = as.numeric(gsub(",", ".", SWEPT_AREA)),
         community_Ic = Ic_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 / 24)) |> 
  group_by(HEX_ID, HAUL_NUMBER, YEAR, MONTH, HAULING_TIME, MEAN_LONGITUDE_DEC, MEAN_LATITUDE_DEC, Trophic_category) |> 
  summarise(community_Ic = sum(community_Ic), .groups = "drop") |> 
  dplyr::filter(Trophic_category %in% c("benthivorous", "generalist", "planktivorous")) |> 
  mutate(Ic_plank = case_when(Trophic_category == "planktivorous" ~ community_Ic,
                              Trophic_category == "generalist" ~ 0.25 * community_Ic, TRUE ~ 0),
         Ic_benthivorous = case_when(Trophic_category == "benthivorous"  ~ community_Ic,
                                     Trophic_category == "generalist" ~ 0.25 * community_Ic, TRUE ~ 0)) |> 
  group_by(HEX_ID, HAUL_NUMBER, YEAR, MONTH, HAULING_TIME, MEAN_LONGITUDE_DEC, MEAN_LATITUDE_DEC) |> 
  summarise(Ic_plank = sum(Ic_plank), Ic_benthivorous = sum(Ic_benthivorous), .groups = "drop") |> 
  mutate(top1_community_Ic_plank = rank(-Ic_plank, ties.method = "first") <= 127,
         top1_community_Ic_benthivorous = rank(-Ic_benthivorous, ties.method = "first") <= 127) 

# Merge datasets at the community and trophic guild levels
Medit_FunCatch_without_NA_community = Medit_FunCatch_without_NA_community |> 
  full_join(Medit_FunCatch_without_NA_trophic) |> 
  relocate(Ic_plank, Ic_benthivorous, .after = community_Gc)
  
# Quick Viz
i = 321
{ID = unique(Medit_FunCatch_without_NA_community$HEX_ID)[i]
(Medit_FunCatch_without_NA_community |> dplyr::filter(HEX_ID == ID) |> 
  ggplot(aes(x = YEAR, y = community_Gc)) + 
    geom_point(size = 4, shape = 21, fill = "orange", color = "black") + 
    ggtitle("Carbon Production via Growth") + theme_classic() +
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18), 
          axis.text  = element_text(size = 16)) +
    scale_x_continuous(name = "", limits = c(1999, 2021), breaks = seq(2000, 2020, 5)) +
    scale_y_continuous(name = expression("Production (gC"~m^-2~d^-1*")"))) +
(Medit_FunCatch_without_NA_community |> dplyr::filter(HEX_ID == ID) |> 
  ggplot(aes(x = YEAR, y = community_Fn)) + 
    geom_point(size = 4, shape = 21, fill = "purple1", color = "black") + 
    ggtitle("Nitrogen Excretion") + theme_classic() +
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18), 
          axis.text  = element_text(size = 16)) +
    scale_x_continuous(name = "", limits = c(1999, 2021), breaks = seq(2000, 2020, 5)) +
    scale_y_continuous(name = expression("N Excretion (gN"~m^-2~d^-1*")")))  +
(Medit_FunCatch_without_NA_community |> dplyr::filter(HEX_ID == ID) |> 
   ggplot(aes(x = YEAR, y = community_Fp)) + 
   geom_point(size = 4, shape = 21, fill = "green4", color = "black") + 
   ggtitle("Phosphorus Excretion") + theme_classic() +
   theme(plot.title = element_text(size = 20),
         axis.title = element_text(size = 18), 
         axis.text  = element_text(size = 16)) +
   scale_x_continuous(name = "", limits = c(1999, 2021), breaks = seq(2000, 2020, 5)) +
   scale_y_continuous(name = expression("P Excretion (gP"~m^-2~d^-1*")")))}

# Spatial analysis Viz
medits_coords <- Medit_FunCatch_without_NA_community |>
  mutate(Longitude = as.numeric(gsub(",", ".", MEAN_LONGITUDE_DEC)),
         Latitude  = as.numeric(gsub(",", ".", MEAN_LATITUDE_DEC))) |> drop_na(HEX_ID) |>
  select(YEAR, Longitude, Latitude, Biomass, top1_Biomass, community_Fn, community_Fp, community_Gc, Ic_plank, 
         Ic_benthivorous, top1_community_Fn, top1_community_Fp, top1_community_Gc, top1_community_Ic_plank,
         top1_community_Ic_benthivorous)

medits_sf <- st_as_sf(medits_coords, coords = c("Longitude", "Latitude"), crs = 4326)
medits_sf_percentile = medits_sf |> mutate(
    Biomass_class       = percentile_class(Biomass),
    community_Fn_class  = percentile_class(community_Fn),
    community_Fp_class  = percentile_class(community_Fp),
    community_Gc_class  = percentile_class(community_Gc),
    community_Icp_class = percentile_class(Ic_plank),
    community_Icb_class = percentile_class(Ic_benthivorous))
land <- ne_countries(scale = "medium", returnclass = "sf")

medits_sf_percentile = medits_sf_percentile |> 
  mutate(Biomass_class = recode(medits_sf_percentile$Biomass_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"),
  community_Gc_class   = recode(medits_sf_percentile$community_Gc_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"),
  community_Fn_class   = recode(medits_sf_percentile$community_Fn_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"),
  community_Fp_class   = recode(medits_sf_percentile$community_Fp_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"),
  community_Icp_class   = recode(medits_sf_percentile$community_Icp_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"),
  community_Icb_class   = recode(medits_sf_percentile$community_Icb_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"))

# Plots
Spatial_Biomass <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, Biomass_class == "Intermediate"),
          aes(fill = Biomass_class, color = Biomass_class, size = Biomass_class, shape = Biomass_class)) +
  geom_sf(data = subset(medits_sf_percentile, Biomass_class %in% c("High", "Low")),
          aes(fill = Biomass_class, color = Biomass_class, size = Biomass_class, shape = Biomass_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#FF968A", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + ggtitle("Fish Biomass") +
  labs(fill  = "Fish Biomass Level", color = "Fish Biomass Level", size  = "Fish Biomass Level", shape = "Fish Biomass Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Spatial_Production <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, community_Gc_class == "Intermediate"),
          aes(fill = community_Gc_class, color = community_Gc_class, size = community_Gc_class, shape = community_Gc_class)) +
  geom_sf(data = subset(medits_sf_percentile, community_Gc_class %in% c("High", "Low")),
          aes(fill = community_Gc_class, color = community_Gc_class, size = community_Gc_class, shape = community_Gc_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#B4CBF0", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + ggtitle("Fish Production") +
  labs(fill = "Fish Production Level", color = "Fish Production Level", size = "Fish Production Level", 
       shape = "Fish Production Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Spatial_Nitrogen <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, community_Fn_class == "Intermediate"),
          aes(fill = community_Fn_class, color = community_Fn_class, size = community_Fn_class, shape = community_Fn_class)) +
  geom_sf(data = subset(medits_sf_percentile, community_Fn_class %in% c("High", "Low")),
          aes(fill = community_Fn_class, color = community_Fn_class, size = community_Fn_class, shape = community_Fn_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#A3B79C", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + ggtitle("Nitrogen Fish Excretion") +
  labs(fill = "Nitrogen Fish Excretion Level", color = "Nitrogen Fish Excretion Level", size = "Nitrogen Fish Excretion Level", 
       shape = "Nitrogen Fish Excretion Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Spatial_Phosphorus <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, community_Fp_class == "Intermediate"),
          aes(fill = community_Fp_class, color = community_Fp_class, size = community_Fp_class, shape = community_Fp_class)) +
  geom_sf(data = subset(medits_sf_percentile, community_Fp_class %in% c("High", "Low")),
          aes(fill = community_Fp_class, color = community_Fp_class, size = community_Fp_class, shape = community_Fp_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#FFF1BA", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + ggtitle("Phosphorus Fish Excretion") +
  labs(fill = "Phosphorus Fish Excretion Level", color = "Phosphorus Fish Excretion Level", size = "Phosphorus Fish Excretion Level", 
       shape = "Phosphorus Fish Excretion Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")
  
Spatial_Planktivory <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, community_Icp_class == "Intermediate"),
          aes(fill = community_Icp_class, color = community_Icp_class, size = community_Icp_class, shape = community_Icp_class)) +
  geom_sf(data = subset(medits_sf_percentile, community_Icp_class %in% c("High", "Low")),
          aes(fill = community_Icp_class, color = community_Icp_class, size = community_Icp_class, shape = community_Icp_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#CCA9DD", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + ggtitle("Fish Planktivory") +
  labs(fill = "Fish planktivory Level", color = "Fish planktivory Level", size = "Fish planktivory Level", 
       shape = "Fish planktivory Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Spatial_Benthivory <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, community_Icb_class == "Intermediate"),
          aes(fill = community_Icb_class, color = community_Icb_class, size = community_Icb_class, shape = community_Icb_class)) +
  geom_sf(data = subset(medits_sf_percentile, community_Icb_class %in% c("High", "Low")),
          aes(fill = community_Icb_class, color = community_Icb_class, size = community_Icb_class, shape = community_Icb_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#FAC898", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + ggtitle("Fish Benthivory") +
  labs(fill = "Fish benthivory Level", color = "Fish benthivory Level", size = "Fish benthivory Level", 
       shape = "Fish benthivory Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_1 = Spatial_Biomass + Spatial_Production + Spatial_Nitrogen + Spatial_Phosphorus +
  Spatial_Planktivory + Spatial_Benthivory + 
  plot_layout(guides = "collect", ncol = 2) & theme(legend.position = "none")

### Add Multifunction index from Schiettekatte et al., 2021
medits_sf_percentile <- medits_sf_percentile |>
  mutate(across(all_of(c("community_Fn", "community_Fp", "community_Gc", "Ic_plank", "Ic_benthivorous")), 
                ~ (.x - min(.x, na.rm = TRUE)) / 
                  (max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)),
                .names = "{.col}_norm")) |> rowwise() |>
  mutate(Multifunctionality = mean(c_across(ends_with("_norm")), na.rm = TRUE)) |> ungroup() |> 
  mutate(Multifunctionality_class = percentile_class(Multifunctionality),
         Multifunctionality_class = recode(Multifunctionality_class, 
  "top_5%" = "Intermediate", "tail_5%" = "Intermediate", "mid" = "Intermediate", "top_1%" = "High", "tail_1%" = "Low"))

Spatial_Multi <- ggplot() +
  geom_sf(data = subset(medits_sf_percentile, Multifunctionality_class == "Intermediate"),
          aes(fill = Multifunctionality_class, color = Multifunctionality_class, size = Multifunctionality_class, 
              shape = Multifunctionality_class)) +
  geom_sf(data = subset(medits_sf_percentile, Multifunctionality_class %in% c("High", "Low")),
          aes(fill = Multifunctionality_class, color = Multifunctionality_class, size = Multifunctionality_class, 
              shape = Multifunctionality_class)) +
  scale_shape_manual(values = c("High" = 21, "Intermediate" = 21, "Low" = 20)) +
  scale_fill_manual(values = c("High" = "#CAE9F5", "Intermediate" = "grey", "Low" = "black")) +
  scale_color_manual(values = c("High" = "black", "Intermediate" = "grey50", "Low" = "black")) +
  scale_size_manual(values = c("High" = 4, "Intermediate" = 3, "Low" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle("Fish Multifunctionality") +
  labs(fill = "Fish multifunctionality Level", color = "Fish multifunctionality Level", 
       size = "Fish multifunctionality Level", shape = "Fish multifunctionality Level") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20, hjust = 0),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_1_tot = ((Spatial_Biomass + Spatial_Production + Spatial_Nitrogen + 
                  Spatial_Phosphorus + Spatial_Planktivory + Spatial_Benthivory) + 
  plot_layout(ncol = 2, guides = "collect")) / plot_spacer() / Spatial_Multi + 
  plot_layout(heights = c(3, 0, 1.68)) & theme(legend.position = "none")

## Let's delineate function distributions
Figure_2A <- ggplot(medits_sf_percentile, aes(x = Biomass)) +
  geom_histogram(binwidth = 0.2, fill = "#FF968A", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("Biomass ("*g~m^{-2}~d^{-1}*")"),
                         y = "Nb of Obs.") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$Biomass), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$Biomass), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#f4cccc") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$Biomass))), y = 0,
             inherit.aes = F, size = 5, fill = "#FF968A", color = "#f4cccc", shape = 21) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2B <- ggplot(medits_sf_percentile, aes(x = community_Fn)) +
  geom_histogram(binwidth = 0.2, fill = "#A3B79C", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("Nitrogen excretion ("*gN~m^{-2}~d^{-1}*")"), y = "") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$community_Fn), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$community_Fn), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#d9ead3") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$community_Fn))), y = 0,
             inherit.aes = F, size = 5, fill = "#A3B79C", color = "#d9ead3", shape = 21) + 
  scale_y_continuous(limits = c(0,3650), breaks = seq(0,3000,1000)) +
  add_phylopic(name = "Mola mola", x = 1e-02, y = 3000, height = 1500)+
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2C <- ggplot(medits_sf_percentile, aes(x = community_Fp)) +
  geom_histogram(binwidth = 0.2, fill = "#FFF1BA", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("Phosphorus excretion ("*gP~m^{-2}~d^{-1}*")"),
                         y = "Nb of Obs.") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$community_Fp), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$community_Fp), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#fff2cc") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$community_Fp))), y = 0,
             inherit.aes = F, size = 5, fill = "#FFF1BA", color = "#fff2cc", shape = 21) +  
  scale_y_continuous(limits = c(0,3650), breaks = seq(0,3000,1000)) +
  add_phylopic(name = "Citharoides macrolepis", x = 2.5e-03, y = 3000, height = 1000) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2D <- ggplot(medits_sf_percentile, aes(x = community_Gc)) +
  geom_histogram(binwidth = 0.2, fill = "#B4CBF0", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("production ("*gC~m^{-2}~d^{-1}*")"), y = "") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$community_Gc), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$community_Gc), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#cfe2f3") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$community_Gc))), y = 0,
             inherit.aes = F, size = 5, fill = "#B4CBF0", color = "#cfe2f3", shape = 21) +  
  scale_y_continuous(limits = c(0,3650), breaks = seq(0,3000,1000)) +
  add_phylopic(name = "Archosargus rhomboidalis", x = 2e-02, y = 3000, height = 1100) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2E <- ggplot(medits_sf_percentile, aes(x = Ic_plank)) +
  geom_histogram(binwidth = 0.2, fill = "#CCA9DD", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("Planktivory ("*gC~m^{-2}~d^{-1}*")"),
                         y = "Nb of Obs.") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$Ic_plank), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$Ic_plank), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#d9d2e9") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$Ic_plank))), y = 0,
             inherit.aes = F, size = 5, fill = "#CCA9DD", color = "#d9d2e9", shape = 21) +  
  scale_y_continuous(limits = c(0,3650), breaks = seq(0,3000,1000)) +
  add_phylopic(name = "Syngnathus typhle", x = 5e-02, y = 3000, height = 600) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2F <- ggplot(medits_sf_percentile, aes(x = Ic_benthivorous)) +
  geom_histogram(binwidth = 0.2, fill = "#FAC898", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("Benthivory ("*gC~m^{-2}~d^{-1}*")"), y = "") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$Ic_benthivorous), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$Ic_benthivorous), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#fff2cc") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$Ic_benthivorous))), y = 0,
             inherit.aes = F, size = 5, fill = "#FAC898", color = "#fff2cc", shape = 21) +
  scale_y_continuous(limits = c(0,3650), breaks = seq(0,3000,1000)) +
  add_phylopic(name = "Raja asterias", x = 9e-02, y = 3000, height = 1300) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2G <- ggplot(medits_sf_percentile, aes(x = Multifunctionality)) +
  geom_histogram(binwidth = 0.2, fill = "#CAE9F5", color = "black") + scale_x_log10() + 
  theme_minimal() + labs(x = expression("Multifunctionality"),
                         y = "Nb of Obs.") +
  geom_linerange(aes(xmin = quantile(na.omit(medits_sf_percentile$Multifunctionality), 0.01),
                     xmax = quantile(na.omit(medits_sf_percentile$Multifunctionality), 0.99), y = 0),
                 inherit.aes = F, size = 2, color = "#cfe2f3") +
  geom_point(aes(x = median(na.omit(medits_sf_percentile$Multifunctionality))), y = 0,
             inherit.aes = F, size = 5, fill = "#CAE9F5", color = "#cfe2f3", shape = 21) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2_tot = ((Figure_2A + Figure_2B + Figure_2C + Figure_2D + Figure_2E + Figure_2F) + 
                  plot_layout(ncol = 2, guides = "collect")) / plot_spacer() / Figure_2G + 
  plot_layout(heights = c(3, 0, 1.68)) & theme(legend.position = "none")

#### Export the data  ----
## Data
# save(Medit_FunCatch_without_NA, file = "Outputs/dat_proc/Medit_FunCatch_without_NA.RData")
# save(medits_sf_percentile, file = "Outputs/dat_proc/medits_sf_percentile.Rdata")
  
## Figures
ggsave(Figure_1_tot, filename = "Figure_1.png", path = "Outputs/", device = "png", width = 12,  height = 14, dpi = 300)  
ggsave(Figure_2_tot, filename = "Figure_2.png", path = "Outputs/", device = "png", width = 10,  height = 12, dpi = 300)  