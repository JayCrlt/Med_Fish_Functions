#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Hexagonal_grid     <- st_read("Data/Grid_0-1000m_Med.shp")
Guilds             <- readxl::read_xlsx("Data/Guilds_MED.xlsx") |> 
  mutate(SPECIES = str_replace_all(SPECIES, "\u00A0", " "))

## Charge from previous scripts
load("Outputs/dat_proc/Med_all.RData")
load("Outputs/dat_proc/Medit_Temp.RData")
load("Outputs/dat_proc/medits_sf_percentile.Rdata")
load("Outputs/dat_proc/Medit_FunCatch_without_NA.RData")

## Color palette and map
land       <- ne_countries(scale = "medium", returnclass = "sf")

## Functions
source("Scripts/00_functions_script.R")

mass = 1000
Med_all_th <- Med_all |> mutate(
  SMR = 0.002 * mass^0.836, MMR = 0.006 * mass^0.779,    
  f0 = SMR / mass^0.836 * exp(-0.138 * (1 / (25 + 273.15) - 1 / (25 + 273.15)) / 8.617e-5),
  f0 = if_else(is.nan(f0), mean(f0[!is.nan(f0)], na.rm = TRUE), f0),
  theta  = (SMR + MMR) / (2 * SMR), TL = (mass / lwa)^(1 / lwb), temp = 25)

## FishFlux loop workflow
# Prepare storing vectors
model         <- list()
param_dataset <- list()
# Prepare the dataset
Med_all_th = Med_all_th |> mutate(Fn_mean = NA, Fn_Q1 = NA, Fn_Q3 = NA, Fp_mean = NA, Fp_Q1 = NA, Fp_Q3 = NA, 
  Gc_mean = NA, Gc_Q1 = NA, Gc_Q3 = NA, Ic_mean = NA, Ic_Q1 = NA, Ic_Q3 = NA)

# Initiate the FishFlux formatting
for(i in 1:length(Med_all_th$Species)){
  param_dataset[[i]] <- list(
    ac_m    = Med_all_th$ac[i],
    ac_sd   = Med_all_th$ac[i] * 0.05,
    an_m    = Med_all_th$an[i],
    an_sd   = Med_all_th$an[i] * 0.05,
    ap_m    = Med_all_th$ap[i],
    ap_sd   = Med_all_th$ap[i] * 0.05,
    r_m     = Med_all_th$r[i],
    f0_m    = Med_all_th$f0[i],
    theta_m = Med_all_th$theta[i],
    Dc_m    = Med_all_th$Dc[i],
    Dc_sd   = Med_all_th$Dc[i] * 0.05,
    Dn_m    = Med_all_th$Dn[i],
    Dn_sd   = Med_all_th$Dn[i] * 0.05,
    Dp_m    = Med_all_th$Dp[i],
    Dp_sd   = Med_all_th$Dp[i] * 0.05,
    k_m     = Med_all_th$K[i],
    linf_m  = Med_all_th$linf[i],
    lwa_m   = Med_all_th$lwa[i],
    lwb_m   = Med_all_th$lwb[i],
    Qc_m    = Med_all_th$Qc[i],
    Qc_sd   = Med_all_th$Qc[i] * 0.05,
    Qn_m    = Med_all_th$Qn[i],
    Qn_sd   = Med_all_th$Qn[i] * 0.05,
    Qp_m    = Med_all_th$Qp[i],
    Qp_sd   = Med_all_th$Qp[i] * 0.05,
    t0_m    = Med_all_th$t0[i],
    h_m     = Med_all_th$h[i],
    F0nz_m  = Med_all_th$F0Nz[i],
    F0pz_m  = Med_all_th$F0Pz[i],
    mdw_m   = Med_all_th$mdw[i],
    v_m     = Med_all_th$temp[i],
    alpha_m = 0.836)}

# Compile data from FishFlux
pb <- txtProgressBar(min = 0, max = length(Med_all_th$Species), style = 3)
for(i in 1: length(Med_all_th$Species)){
  tryCatch({
    invisible(capture.output({
      model[[i]] <- extract(cnp_model_mcmc(TL = Med_all_th$TL[i], 
                                           param = param_dataset[[i]]), c("Fn","Fp", "Gc", "Ic"))}))
    Med_all_th$Fn_mean[i] <- model[[i]]$Fn_mean
    Med_all_th$Fn_Q1[i]   <- model[[i]]$`Fn_2.5%`
    Med_all_th$Fn_Q3[i]   <- model[[i]]$`Fn_97.5%`
    Med_all_th$Fp_mean[i] <- model[[i]]$Fp_mean
    Med_all_th$Fp_Q1[i]   <- model[[i]]$`Fp_2.5%`
    Med_all_th$Fp_Q3[i]   <- model[[i]]$`Fp_97.5%`
    Med_all_th$Gc_mean[i] <- model[[i]]$Gc_mean
    Med_all_th$Gc_Q1[i]   <- model[[i]]$`Gc_2.5%`
    Med_all_th$Gc_Q3[i]   <- model[[i]]$`Gc_97.5%`
    Med_all_th$Ic_mean[i] <- model[[i]]$Ic_mean
    Med_all_th$Ic_Q1[i]   <- model[[i]]$`Ic_2.5%`
    Med_all_th$Ic_Q3[i]   <- model[[i]]$`Ic_97.5%`
  }, error = function(e){
    message("Error in row ", i)
    Med_all_th$Fn_mean[i] <- NA
    Med_all_th$Fn_Q1[i]   <- NA
    Med_all_th$Fn_Q3[i]   <- NA
    Med_all_th$Fp_mean[i] <- NA
    Med_all_th$Fp_Q1[i]   <- NA
    Med_all_th$Fp_Q3[i]   <- NA
    Med_all_th$Gc_mean[i] <- NA
    Med_all_th$Gc_Q1[i]   <- NA
    Med_all_th$Gc_Q3[i]   <- NA
    Med_all_th$Ic_mean[i] <- NA
    Med_all_th$Ic_Q1[i]   <- NA
    Med_all_th$Ic_Q3[i]   <- NA
  })
  setTxtProgressBar(pb, i)}

Med_all_th_sum <- Med_all_th |> group_by(Family) |> summarise(Fp_family_mean = mean(Fp_mean),
                                                              Fp_family_sd   = sd(Fp_mean),
                                                              Fn_family_mean = mean(Fn_mean),
                                                              Fn_family_sd   = sd(Fn_mean),
                                                              Gc_family_mean = mean(Gc_mean),
                                                              Gc_family_sd   = sd(Gc_mean),
                                                              Ic_family_mean = mean(Ic_mean),
                                                              Ic_family_sd   = sd(Ic_mean)) |> mutate(uuid = NA)

# Look for icons using Rphylopic
for (i in 1:length(Med_all_th_sum$uuid)) {
  Med_all_th_sum$uuid[i] <- tryCatch(rphylopic::get_uuid(name = gsub("_", " ", Med_all_th_sum$Family[i]), n = 1),
                                     error = function(e) NA_character_)}

# Complete the dataset manually for famiy not found
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Argentinidae"]      = "e9f4e72d-377e-44ed-8722-0362fbd589ad"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Aulopidae"]         = "f2396149-3019-45a6-b87c-9ac5bfeba8b1"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Bothidae"]          = "2f31392b-b45b-4971-8e70-1592b9554da8"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Bramidae"]          = "0824b100-e616-4a6c-983e-4d8d0f0282ba"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Bythitidae"]        = "7a0c849a-9ca0-4411-abbb-357a43ed9944"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Gonostomatidae"]    = "5400612d-8088-470b-a7b0-36ecde18bbaf"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Anthiadidae"]       = "eb6f54bf-3182-4d3a-9f16-602759c76db0"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Caproidae"]         = "922483bb-ae6a-43e2-b3ff-45dc54ebbd87"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Synaphobranchidae"] = "33d68cc1-4702-4d9e-8e23-5aa07acaf3ee"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Nettastomatidae"]   = "f2f4c62b-e359-45c5-ad69-29041e5f2ceb"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Scomberesocidae"]   = "3b3f7a2f-8997-4969-956c-201de271707b"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Chlorophthalmidae"] = "1ac790b5-d34f-4568-9f20-e83060b48abb"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Cepolidae"]         = "23c0bb71-c7f1-4ce2-9bf9-e1c8947f2501"
Med_all_th_sum$uuid[Med_all_th_sum$Family == "Liparidae"]         = "380c8d99-e483-4473-a6aa-fa3668e58279"


Phosphorus  = Med_all_th |> left_join(Med_all_th_sum) |> distinct(Family, Fp_family_mean) |> 
  arrange(desc(Fp_family_mean)) |> slice_head(n = 10) |> 
  left_join(Med_all_th_sum |> dplyr::select(Family, uuid)) |> 
  rename(Mean_rate = Fp_family_mean) |> mutate(variable = "Fp")
Nitrogen    = Med_all_th |> left_join(Med_all_th_sum) |> distinct(Family, Fn_family_mean) |> 
  arrange(desc(Fn_family_mean)) |> slice_head(n = 10) |> 
  left_join(Med_all_th_sum |> dplyr::select(Family, uuid)) |> 
  rename(Mean_rate = Fn_family_mean) |> mutate(variable = "Fn")
Carbon      = Med_all_th |> left_join(Med_all_th_sum) |> distinct(Family, Gc_family_mean) |> 
  arrange(desc(Gc_family_mean)) |> slice_head(n = 10) |> 
  left_join(Med_all_th_sum |> dplyr::select(Family, uuid)) |> 
  rename(Mean_rate = Gc_family_mean) |> mutate(variable = "Gc")
Planktivory = Med_all_th |> left_join(Med_all_th_sum) |> 
  left_join(Guilds |> mutate(Species = gsub(" ", "_", SPECIES)) |> dplyr::select(Species, Trophic_category)) |> 
  dplyr::filter(Trophic_category %in% c("planktivorous", "generalist")) |> 
  dplyr::mutate(Ic_family_mean = dplyr::case_when(Trophic_category == "generalist" ~ Ic_family_mean / 4,
                                           Trophic_category == "planktivorous" ~ Ic_family_mean, 
                                           TRUE ~ Ic_family_mean)) |> 
  distinct(Family, Ic_family_mean) |> arrange(desc(Ic_family_mean)) |> slice_head(n = 10) |> 
  left_join(Med_all_th_sum |> dplyr::select(Family, uuid)) |> 
  rename(Mean_rate = Ic_family_mean) |> mutate(variable = "IcP")
Benthivory  = Med_all_th |> left_join(Med_all_th_sum) |> 
  left_join(Guilds |> mutate(Species = gsub(" ", "_", SPECIES)) |> dplyr::select(Species, Trophic_category)) |> 
  dplyr::filter(Trophic_category %in% c("benthivorous", "generalist")) |> 
  dplyr::mutate(Ic_family_mean = dplyr::case_when(Trophic_category == "generalist" ~ Ic_family_mean / 4,
                                                  Trophic_category == "planktivorous" ~ Ic_family_mean, 
                                                  TRUE ~ Ic_family_mean)) |> 
  distinct(Family, Ic_family_mean) |> arrange(desc(Ic_family_mean)) |> slice_head(n = 10) |> 
  left_join(Med_all_th_sum |> dplyr::select(Family, uuid)) |> 
  rename(Mean_rate = Ic_family_mean) |> mutate(variable = "IcB")

Family_rates <- rbind(Phosphorus, Nitrogen, Carbon, Planktivory, Benthivory) |> 
  mutate(ID = rep(seq(1,10,1), 5))


Figure_2A = Family_rates |> 
  filter(variable == "Fn") |> 
  ggplot(aes(x = reorder(ID, Mean_rate), y = Mean_rate)) +
  geom_segment(aes(xend = reorder(ID, Mean_rate), y = 0, yend = Mean_rate), color = "grey60", size = 1.2) +
  geom_point(size = 5, fill = "#A3B79C", color = "black", shape = 21) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_flip() + theme_classic() +  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(x = NULL, y = expression(atop("Nitrogen excretion", "(gN kg"^{-1}~d^{-1}*")"))) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2B = Family_rates |> 
  filter(variable == "Fp") |> 
  ggplot(aes(x = reorder(ID, Mean_rate), y = Mean_rate)) +
  geom_segment(aes(xend = reorder(ID, Mean_rate), y = 0, yend = Mean_rate), color = "grey60", size = 1.2) +
  geom_point(size = 5, fill = "#FFF1BA", color = "black", shape = 21) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_flip() + theme_classic() +  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(x = NULL, y = expression(atop("Phosphorus excretion", "(gP kg"^{-1}~d^{-1}*")"))) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2C = Family_rates |> 
  filter(variable == "Gc") |> 
  ggplot(aes(x = reorder(ID, Mean_rate), y = Mean_rate)) +
  geom_segment(aes(xend = reorder(ID, Mean_rate), y = 0, yend = Mean_rate), color = "grey60", size = 1.2) +
  geom_point(size = 5, fill = "#B4CBF0", color = "black", shape = 21) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_flip() + theme_classic() +  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(x = NULL, y = expression(atop("Production", "(gC kg"^{-1}~d^{-1}*")"))) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2D = Family_rates |> 
  filter(variable == "IcP") |> 
  ggplot(aes(x = reorder(ID, Mean_rate), y = Mean_rate)) +
  geom_segment(aes(xend = reorder(ID, Mean_rate), y = 0, yend = Mean_rate), color = "grey60", size = 1.2) +
  geom_point(size = 5, fill = "#CCA9DD", color = "black", shape = 21) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_flip() + theme_classic() +  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(x = NULL, y = expression(atop("Planktivory", "(gC kg"^{-1}~d^{-1}*")"))) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2E = Family_rates |> 
  filter(variable == "IcB") |> 
  ggplot(aes(x = reorder(ID, Mean_rate), y = Mean_rate)) +
  geom_segment(aes(xend = reorder(ID, Mean_rate), y = 0, yend = Mean_rate), color = "grey60", size = 1.2) +
  geom_point(size = 5, fill = "#FAC898", color = "black", shape = 21) +
  scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  coord_flip() + theme_classic() +  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(x = NULL, y = expression(atop("Benthivory", "(gC kg"^{-1}~d^{-1}*")"))) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_2_tot = (Figure_2A + Figure_2B + Figure_2C + Figure_2D + Figure_2E) +
                  plot_layout(ncol = 3)

#### Export the data  ----
## Figures
ggsave(Figure_2_tot, filename = "Figure_2.png", path = "Outputs/", device = "png", width = 10,  height = 8, dpi = 300)  