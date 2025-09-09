#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Hexagonal_grid     <- st_read("Data/Grid 0-1000m_WMED/Grid 0-1000m_WMED.shp")

## Charge from previous scripts
load("Outputs/dat_proc/Western_Med.RData")
load("Outputs/dat_proc/Medit_Temp.RData")

## Functions
source("Scripts/00_functions_script.R")

# add spatial hex info
medits_sf               <- st_as_sf(data.frame(Longitude = as.numeric(gsub(",", ".", MEDIT_TEMP$MEAN_LONGITUDE_DEC)), 
                                               Latitude = as.numeric(gsub(",", ".", MEDIT_TEMP$MEAN_LATITUDE_DEC))), 
                                    coords = c("Longitude", "Latitude"), crs = 4326)
medits_sf_hex           <- st_join(medits_sf, Hexagonal_grid["grid.id"], left = TRUE)
MEDIT_TEMP$HEX_ID       <- medits_sf_hex$grid.id
MEDIT_TEMP$MASS_IND     <- MEDIT_TEMP$TOTAL_WEIGHT_IN_THE_HAUL / MEDIT_TEMP$TOTAL_NUMBER_IN_THE_HAUL

# Clean data according to biomass vs species
Medits_total            <- MEDIT_TEMP |> dplyr::filter(TOTAL_WEIGHT_IN_THE_HAUL != 0 & TOTAL_NUMBER_IN_THE_HAUL != 0) |> 
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
Medit_Western_FunCatch = Medit_Western_FunCatch |> 
  mutate(Fn_mean = NA, Fn_Q1 = NA, Fn_Q3 = NA, Fp_mean = NA, Fp_Q1 = NA, Fp_Q3 = NA, 
         Gc_mean = NA, Gc_Q1 = NA, Gc_Q3 = NA, Ic_mean = NA, Ic_Q1 = NA, Ic_Q3 = NA)

# Initiate the FishFlux formatting
for(i in 1:length(Medit_Western_FunCatch$HEX_ID)){
  param_dataset[[i]] <- list(
    ac_m    = Medit_Western_FunCatch$ac[i],
    ac_sd   = Medit_Western_FunCatch$ac[i] * 0.05,
    an_m    = Medit_Western_FunCatch$an[i],
    an_sd   = Medit_Western_FunCatch$an[i] * 0.05,
    ap_m    = Medit_Western_FunCatch$ap[i],
    ap_sd   = Medit_Western_FunCatch$ap[i] * 0.05,
    r_m     = Medit_Western_FunCatch$r[i],
    f0_m    = Medit_Western_FunCatch$f0[i],
    theta_m = Medit_Western_FunCatch$theta[i],
    Dc_m    = Medit_Western_FunCatch$Dc[i],
    Dc_sd   = Medit_Western_FunCatch$Dc[i] * 0.05,
    Dn_m    = Medit_Western_FunCatch$Dn[i],
    Dn_sd   = Medit_Western_FunCatch$Dn[i] * 0.05,
    Dp_m    = Medit_Western_FunCatch$Dp[i],
    Dp_sd   = Medit_Western_FunCatch$Dp[i] * 0.05,
    k_m     = Medit_Western_FunCatch$K[i],
    linf_m  = Medit_Western_FunCatch$linf[i],
    lwa_m   = Medit_Western_FunCatch$lwa[i],
    lwb_m   = Medit_Western_FunCatch$lwb[i],
    Qc_m    = Medit_Western_FunCatch$Qc[i],
    Qc_sd   = Medit_Western_FunCatch$Qc[i] * 0.05,
    Qn_m    = Medit_Western_FunCatch$Qn[i],
    Qn_sd   = Medit_Western_FunCatch$Qn[i] * 0.05,
    Qp_m    = Medit_Western_FunCatch$Qp[i],
    Qp_sd   = Medit_Western_FunCatch$Qp[i] * 0.05,
    t0_m    = Medit_Western_FunCatch$t0[i],
    h_m     = Medit_Western_FunCatch$h[i],
    F0nz_m  = Medit_Western_FunCatch$F0Nz[i],
    F0pz_m  = Medit_Western_FunCatch$F0Pz[i],
    mdw_m   = Medit_Western_FunCatch$mdw[i],
    v_m     = Medit_Western_FunCatch$TEMP[i],
    alpha_m = 0.836)}

# Compile data from FishFlux
pb <- txtProgressBar(min = 0, max = length(Medit_Western_FunCatch$HEX_ID), style = 3)
for(i in 1: length(Medit_Western_FunCatch$HEX_ID)){
  tryCatch({
    invisible(capture.output({
  model[[i]] <- extract(cnp_model_mcmc(TL = Medit_Western_FunCatch$TL[i], param = param_dataset[[i]]), c("Fn","Fp", "Gc", "Ic"))}))
  Medit_Western_FunCatch$Fn_mean[i] <- model[[i]]$Fn_mean
  Medit_Western_FunCatch$Fn_Q1[i]   <- model[[i]]$`Fn_2.5%`
  Medit_Western_FunCatch$Fn_Q3[i]   <- model[[i]]$`Fn_97.5%`
  Medit_Western_FunCatch$Fp_mean[i] <- model[[i]]$Fp_mean
  Medit_Western_FunCatch$Fp_Q1[i]   <- model[[i]]$`Fp_2.5%`
  Medit_Western_FunCatch$Fp_Q3[i]   <- model[[i]]$`Fp_97.5%`
  Medit_Western_FunCatch$Gc_mean[i] <- model[[i]]$Gc_mean
  Medit_Western_FunCatch$Gc_Q1[i]   <- model[[i]]$`Gc_2.5%`
  Medit_Western_FunCatch$Gc_Q3[i]   <- model[[i]]$`Gc_97.5%`
  Medit_Western_FunCatch$Ic_mean[i] <- model[[i]]$Ic_mean
  Medit_Western_FunCatch$Ic_Q1[i]   <- model[[i]]$`Ic_2.5%`
  Medit_Western_FunCatch$Ic_Q3[i]   <- model[[i]]$`Ic_97.5%`
  }, error = function(e){
    message("Error in row ", i)
    Medit_Western_FunCatch$Fn_mean[i] <- NA
    Medit_Western_FunCatch$Fn_Q1[i]   <- NA
    Medit_Western_FunCatch$Fn_Q3[i]   <- NA
    Medit_Western_FunCatch$Fp_mean[i] <- NA
    Medit_Western_FunCatch$Fp_Q1[i]   <- NA
    Medit_Western_FunCatch$Fp_Q3[i]   <- NA
    Medit_Western_FunCatch$Gc_mean[i] <- NA
    Medit_Western_FunCatch$Gc_Q1[i]   <- NA
    Medit_Western_FunCatch$Gc_Q3[i]   <- NA
    Medit_Western_FunCatch$Ic_mean[i] <- NA
    Medit_Western_FunCatch$Ic_Q1[i]   <- NA
    Medit_Western_FunCatch$Ic_Q3[i]   <- NA
  })
  setTxtProgressBar(pb, i)}

# There was significant numbers of errors due to missing observed temperature, it has been fixed thanks to R_script_02
Medit_Western_FunCatch_without_NA = Medit_Western_FunCatch |> 
  anti_join(Medit_Western_FunCatch |> group_by(HEX_ID, HAUL_NUMBER, YEAR, MONTH, HAULING_TIME) |> 
              summarise(NaN_count = sum(is.nan(TEMP)), Total_count = n(), .groups = "drop") |> 
              filter(NaN_count > 0) |> select(HEX_ID, YEAR), by = c("HEX_ID", "YEAR")) # Remove 257315 - 257243 = -72 obs (< -0.03%) 

# Community computation
Medit_Western_FunCatch_without_NA_community = Medit_Western_FunCatch_without_NA |> 
  mutate(SWEPT_AREA = as.numeric(gsub(",", ".", SWEPT_AREA)),
         Biomass = TOTAL_WEIGHT_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 * 24),
         community_Fn = Fn_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 * 24),
         community_Fp = Fp_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 * 24),
         community_Gc = Gc_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 * 24),
         community_Ic = Ic_mean * TOTAL_NUMBER_IN_THE_HAUL / HAUL_AREA / (HAUL_DURATION / 60 * 24)) |> 
  group_by(HEX_ID, HAUL_NUMBER, YEAR, MONTH, HAULING_TIME, MEAN_LONGITUDE_DEC, MEAN_LATITUDE_DEC) |> 
  summarise(Biomass = sum(Biomass),
            community_Fn = sum(community_Fn),
            community_Fp = sum(community_Fp),
            community_Gc = sum(community_Gc),
            community_Ic = sum(community_Ic),
            .groups = "drop") |> 
  mutate(top1_Biomass      = rank(-Biomass, ties.method = "first") <= 85,
         top1_community_Fn = rank(-community_Fn, ties.method = "first") <= 85,
         top1_community_Fp = rank(-community_Fp, ties.method = "first") <= 85,
         top1_community_Gc = rank(-community_Gc, ties.method = "first") <= 85)

# Quick Viz
i = 321
{ID = unique(Medit_Western_FunCatch_without_NA_community$HEX_ID)[i]
(Medit_Western_FunCatch_without_NA_community |> dplyr::filter(HEX_ID == ID) |> 
  ggplot(aes(x = YEAR, y = community_Gc)) + 
    geom_point(size = 4, shape = 21, fill = "orange", color = "black") + 
    ggtitle("Carbon Production via Growth") + theme_classic() +
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18), 
          axis.text  = element_text(size = 16)) +
    scale_x_continuous(name = "", limits = c(1999, 2021), breaks = seq(2000, 2020, 5)) +
    scale_y_continuous(name = expression("Production (gC"~m^-2~d^-1*")"))) +
(Medit_Western_FunCatch_without_NA_community |> dplyr::filter(HEX_ID == ID) |> 
  ggplot(aes(x = YEAR, y = community_Fn)) + 
    geom_point(size = 4, shape = 21, fill = "purple1", color = "black") + 
    ggtitle("Nitrogen Excretion") + theme_classic() +
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 18), 
          axis.text  = element_text(size = 16)) +
    scale_x_continuous(name = "", limits = c(1999, 2021), breaks = seq(2000, 2020, 5)) +
    scale_y_continuous(name = expression("N Excretion (gN"~m^-2~d^-1*")")))  +
(Medit_Western_FunCatch_without_NA_community |> dplyr::filter(HEX_ID == ID) |> 
   ggplot(aes(x = YEAR, y = community_Fp)) + 
   geom_point(size = 4, shape = 21, fill = "green4", color = "black") + 
   ggtitle("Phosphorus Excretion") + theme_classic() +
   theme(plot.title = element_text(size = 20),
         axis.title = element_text(size = 18), 
         axis.text  = element_text(size = 16)) +
   scale_x_continuous(name = "", limits = c(1999, 2021), breaks = seq(2000, 2020, 5)) +
   scale_y_continuous(name = expression("P Excretion (gP"~m^-2~d^-1*")")))}

# Spatial analysis Viz
medits_coords <- Medit_Western_FunCatch_without_NA_community |>
  mutate(Longitude = as.numeric(gsub(",", ".", MEAN_LONGITUDE_DEC)),
         Latitude  = as.numeric(gsub(",", ".", MEAN_LATITUDE_DEC))) |>
  select(Longitude, Latitude, Biomass, top1_Biomass, community_Fn, community_Fp, community_Gc,
         top1_community_Fn, top1_community_Fp, top1_community_Gc)
medits_sf <- st_as_sf(medits_coords, coords = c("Longitude", "Latitude"), crs = 4326)
hex_grid <- st_transform(Hexagonal_grid, crs = st_crs(medits_sf))
Hexagonal_grid$Total_Biomass <- sapply(st_intersects(Hexagonal_grid, medits_sf), function(idx) sum(medits_sf$Biomass[idx]))
Hexagonal_grid$Top1_Biomass <- sapply(st_intersects(Hexagonal_grid, medits_sf |>
                                                      filter(top1_Biomass)), length) > 0

for(var in c("community_Fn", "community_Fp", "community_Gc")) {
  Hexagonal_grid[[var]] <- sapply(st_intersects(Hexagonal_grid, medits_sf), function(idx) sum(medits_sf[[var]][idx]))
  Hexagonal_grid[[paste0("top1_", var)]] <- sapply(st_intersects(Hexagonal_grid, medits_sf |>
                                                                   filter(get(paste0("top1_", var)))), length) > 0}
land <- ne_countries(scale = "medium", returnclass = "sf")
hex_marine <- st_difference(Hexagonal_grid, st_union(land)) |>
  arrange(desc(Total_Biomass)) |> distinct(geometry, .keep_all = TRUE)

# Plot
Spatial_Biomass = ggplot() +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_sf(data = hex_marine, aes(fill = Total_Biomass), color = "black", show.legend = FALSE) +
  geom_sf(data = hex_marine |> filter(Top1_Biomass), fill = "red", color = "purple", size = 0.8, show.legend = FALSE) +
  scale_fill_gradient(low = "moccasin", high = "darkorange") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(fill = "Hex Total Biomass", title = "Biomass") +
  coord_sf(xlim = c(-6, 16), ylim = c(35, 45)) +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 18), 
        axis.text  = element_text(size = 16))

Spatial_Production = ggplot() +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_sf(data = hex_marine, aes(fill = community_Gc), color = "black", show.legend = FALSE) +
  geom_sf(data = hex_marine |> filter(top1_community_Gc), fill = "red", color = "purple", size = 0.8, show.legend = FALSE) +
  scale_fill_gradient(low = "moccasin", high = "darkorange") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(fill = "Hex Total Biomass", title = "Carbon Production via Growth") +
  coord_sf(xlim = c(-6, 16), ylim = c(35, 45)) +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 18), 
        axis.text  = element_text(size = 16))

Spatial_Nitrogen = ggplot() +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_sf(data = hex_marine, aes(fill = community_Fn), color = "black", show.legend = FALSE) +
  geom_sf(data = hex_marine |> filter(top1_community_Fn), fill = "red", color = "purple", size = 0.8, show.legend = FALSE) +
  scale_fill_gradient(low = "moccasin", high = "darkorange") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Nitrogen Excretion") +
  coord_sf(xlim = c(-6, 16), ylim = c(35, 45)) +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 18), 
        axis.text  = element_text(size = 16))

Spatial_Phosphorus = ggplot() +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_sf(data = hex_marine, aes(fill = community_Fp), color = "black", show.legend = FALSE) +
  geom_sf(data = hex_marine |> filter(top1_community_Fp), fill = "red", color = "purple", size = 0.8, show.legend = FALSE) +
  scale_fill_gradient(low = "moccasin", high = "darkorange") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(title = "Phosphorus Excretion") +
  coord_sf(xlim = c(-6, 16), ylim = c(35, 45)) +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 18), 
        axis.text  = element_text(size = 16))

# regression for each cell
Figure_1 = Spatial_Biomass + Spatial_Production + Spatial_Nitrogen + Spatial_Phosphorus

#### Export the data  ----
## Data
# save(Medit_Western_FunCatch_without_NA, file = "Outputs/dat_proc/Medit_Western_FunCatch_without_NA.RData")

## Figures
ggsave(Figure_1, filename = "Figure_1.png", path = "Outputs/", device = "png", width = 20,  height = 12, dpi = 300)  