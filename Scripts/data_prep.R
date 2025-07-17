#### Setting up          ----
rm(list = ls())
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata")

## Download data
# Raw data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Medits_total       <- read.delim("Data/TATB_WMED_1999-2021_clean.csv", sep = ";")
Nutrients_FishBase <- readxl::read_excel("Data/Nutrients_FishBase.xlsx") #87 values for C, 206 for N, 46 for P
Metabolic_Rosen    <- read.delim("Data/Rosen_2025.csv", sep = ";")
Metabolic_Nina     <- read.delim("Data/Schiettekatte_2021_Metabolism.csv", sep = ",")
Hexagonal_grid     <- st_read("Data/Grid 0-1000m_WMED/Grid 0-1000m_WMED.shp")

# Nina data CNP diet
cnp_diet           <- read.delim("Data/cnp_diet.csv", sep = ",") |> 
  rename(Species = species) |> 
  left_join(rfishbase::load_taxa(), by = "Species") |> 
  left_join(rfishbase::ecology(), by = "SpecCode") 

## Functions
# Imputational function
impute_trait_phylopars <- function(data, trait_col, phy) {
  trait_vector <- data |>
    filter(!is.na(.data[[trait_col]])) |>
    transmute(Species = gsub(" ", "_", Species), Trait = .data[[trait_col]]) |>
    tibble::deframe()
  trait_data <- data.frame(Species = names(trait_vector), Trait = trait_vector)
  trait_data$Species <- gsub(" ", "_", trait_data$Species)
  matched_species <- intersect(trait_data$Species, phy$tip.label)
  trait_data_clean <- trait_data |> filter(Species %in% matched_species)
  fit <- phylopars(trait_data = trait_data_clean |> rename(species = Species), tree = phy)
  imputed_matrix <- fit$anc_recon[phy$tip.label, , drop = FALSE]
  trait_name <- gsub("[^a-zA-Z]", "", trait_col)  
  estimates <- tibble(Species = rownames(imputed_matrix), Estimated = imputed_matrix[, "Trait"])
  result <- data |> 
    left_join(estimates, by = "Species") |>
    mutate(
      !!paste0(trait_name, "_final") := ifelse(is.na(.data[[trait_col]]), Estimated, .data[[trait_col]]),
      !!paste0(trait_name, "_type") := ifelse(is.na(.data[[trait_col]]), "imputed", "measured")) |>
    select(-Estimated, -all_of(trait_col)) |> rename(!!trait_name := paste0(trait_name, "_final"))
  return(result)}
# Sampling
sample_once <- function(df, n_sample) {
  df |> group_by(YEAR) |> slice_sample(n = n_sample) |>
    mutate(HAUL_DURATION = as.numeric(HAUL_DURATION), DISTANCE = as.numeric(DISTANCE), TOTAL_WEIGHT_DIST_KG_KM_H = 
             (TOTAL_WEIGHT_IN_THE_HAUL / 1000) / (DISTANCE / 1000) / (HAUL_DURATION * 60)) |>
    group_by(YEAR) |> summarise(weight_std = sum(TOTAL_WEIGHT_DIST_KG_KM_H, na.rm = T)) |> ungroup()}
# Function to fill NAs hierarchically for one trait
fill_trait_hierarchy <- function(df, trait) {
  trait_sym <- sym(trait)
  df <- df %>%
    group_by(Genus) %>%
    mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE),!!trait_sym)) %>%
    ungroup() %>% group_by(Family) %>%
    mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE),!!trait_sym)) %>%
    ungroup() %>% mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE), !!trait_sym))
  return(df)}

#### 0. Exploration      ----
# Merge Species_code with Scientific name
Medits_total = Medits_total |> mutate(spp.code = paste(GENUS, SPECIES, sep = "")) |> 
  left_join(sp_code_list, by = "spp.code") |> 
  filter(!str_detect(sci.name, "^NO\\s*")) |> 
  filter(catfau %in% c("Aa", "Ae", "Ao")) # Work only with fish first

# How many species in the Med Sea to play with?
length(unique(Medits_total$sci.name)) # 361 fish out of 1246 species in total

# So let's compile first the data we need to evaluate the gaps using rfishbase
B_Useful_species_list = Medits_total |> group_by(sci.name) |> distinct(sci.name) |> data.frame() |> 
  rename(Species = sci.name) |> 
  left_join(Nutrients_FishBase, by = "Species") |> 
  left_join(rfishbase::load_taxa(), by = "Species") |> 
  left_join(rfishbase::popqb(), by = "SpecCode") |> 
  dplyr::select(Species, Genus, Family.y, `C/DM`, `N/DM`, `P/DM`, PopQB, FoodType) |> 
  rename(Family = Family.y)

# Look at global data
medits_coords           <- data.frame(Longitude = as.numeric(gsub(",", ".", Medits_total$MEAN_LONGITUDE_DEC)), 
                                      Latitude = as.numeric(gsub(",", ".", Medits_total$MEAN_LATITUDE_DEC)))
medits_sf               <- st_as_sf(medits_coords, coords = c("Longitude", "Latitude"), crs = 4326)
hex_grid                <- st_transform(Hexagonal_grid, crs = st_crs(medits_sf))
Hexagonal_grid$ObsCount <- lengths(st_intersects(Hexagonal_grid, medits_sf))
suppressWarnings({land  <- ne_countries(scale = "medium", returnclass = "sf") |> 
  st_crop(xmin = -6, xmax = 36, ymin = 30, ymax = 46)
  hex_marine <- st_difference(Hexagonal_grid |> filter(ObsCount > 0), st_union(land))})
Figure_1 <- leaflet() |> addProviderTiles(providers$Esri.WorldImagery) |>
  addPolygons(data = hex_marine,
              fillColor = ~colorNumeric(palette = "YlOrRd", domain = hex_marine$ObsCount)(ObsCount),
              fillOpacity = 1, color = "red", weight = 1, smoothFactor = 0.2) |>
  addScaleBar(position = "bottomleft") |> addFullscreenControl()

# Catch surveys
Catches_evolution <- map_dfr(1:49, ~ sample_once(Medits_total, 5000) |> mutate(iteration = .x)) |> 
    group_by(YEAR) |> summarise(mean_weight = mean(weight_std, na.rm = T), sd_weight = sd(weight_std, na.rm = T)) |> 
    ggplot(aes(x = YEAR, y = mean_weight)) +
    geom_linerange(aes(ymin = mean_weight - sd_weight, ymax = mean_weight + sd_weight)) + 
    geom_line(linetype = "dotted") + geom_point(shape = 21, fill = "white", size = 3) + 
    theme_classic() + labs(y = expression(atop(paste("Standardized hauling activity (kg.km"^-1, ".h"^-1, ") ± SD"),
                                               "(5000 obs.yr"^-1*" and 50 iterations)")), x = "") +
    theme(axis.title = element_text(size = 16), axis.text  = element_text(size = 14))

#### 1. CNP Body mass    ----
# Distinct species in Nutrients database
Nutrients_FishBase = Nutrients_FishBase |> group_by(Species, Family) |> 
  summarise(`C/DM` = mean(`C/DM`), `N/DM` = mean(`N/DM`), `P/DM` = mean(`P/DM`))
# Standardize species names once
Nutrients_FishBase$Species <- gsub(" ", "_", Nutrients_FishBase$Species)
# Load phylogeny once
phy <- fishtree_phylogeny(species = Nutrients_FishBase$Species)
phy$tip.label <- gsub(" ", "_", phy$tip.label)
nutrients_imputed <- Nutrients_FishBase |>
  impute_trait_phylopars(trait_col = "C/DM", phy = phy) |>
  impute_trait_phylopars(trait_col = "N/DM", phy = phy) |>
  impute_trait_phylopars(trait_col = "P/DM", phy = phy)

C_body = nutrients_imputed |> 
  ggplot() + geom_density_ridges(aes(y = CDM_type, x = CDM*100, fill = CDM_type), 
                                 alpha = .6, bandwidth = 3, show.legend = F) +
  scale_x_continuous(name = "Body mass composition (%C)") +
  scale_y_discrete(name = "Type of data") + theme_classic() +
  scale_fill_manual(values = c("brown", "darkcyan")) +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

N_body = nutrients_imputed |> 
  ggplot() + geom_density_ridges(aes(y = NDM_type, x = NDM*100, fill = NDM_type), 
                                 alpha = .6, bandwidth = .5, show.legend = F) +
  scale_x_continuous(name = "Body mass composition (%N)") +
  scale_y_discrete(name = "Type of data") + theme_classic() +
  scale_fill_manual(values = c("brown", "darkcyan")) +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

P_body = nutrients_imputed |> 
  ggplot() + geom_density_ridges(aes(y = PDM_type, x = PDM*100, fill = PDM_type), 
                                 alpha = .6, bandwidth = .5, show.legend = F) +
  scale_x_continuous(name = "Body mass composition (%P)") +
  scale_y_discrete(name = "Type of data") + theme_classic() +
  scale_fill_manual(values = c("brown", "darkcyan")) +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

Figure_S1 = C_body / N_body / P_body + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

#### 2. CNP diets        ----
data_summary <- cnp_diet |> group_by(Species) |> summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) |> 
  filter(!is.na(FoodTroph), !is.na(c))
model_C      <- lm(c ~ I(FoodTroph - 1) + 0, data = data_summary)
data_summary <- data_summary |> mutate(Fitted = predict(model_C),
    Above = ifelse(c > Fitted, "above", "below"), Color = ifelse(Above == "above", "dodgerblue2", "coral2"))
r2           <- summary(model_C)$r.squared
r2_label     <- paste0("italic(R)^2 == ", format(round(r2, 2), nsmall = 2))
C_diet <- ggplot(data_summary, aes(x = FoodTroph, y = c)) +
  geom_segment(aes(xend = FoodTroph, yend = Fitted), linetype = "dotted") +
  geom_point(aes(color = Above, fill = Color), shape = 21, size = 3, stroke = 1, show.legend = F) +
  scale_fill_identity() +
  scale_color_manual(values = c("above" = "black", "below" = "black")) +
  stat_smooth(method = "lm", formula = y ~ I(x - 1) + 0, color = "mediumpurple", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = r2_label, parse = TRUE,
           hjust = 1.1, vjust = 1.5, size = 5) + theme_classic() +
  scale_x_continuous(name = "Trophic level", limits = c(2, 4.5)) +
  scale_y_continuous(name = "C diet (%)") +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

model_N      <- lm(n ~ I(FoodTroph - 1) + 0, data = data_summary)
data_summary <- data_summary |> mutate(Fitted = predict(model_N),
    Above = ifelse(n > Fitted, "above", "below"), Color = ifelse(Above == "above", "dodgerblue2", "coral2"))
r2           <- summary(model_N)$r.squared
r2_label     <- paste0("italic(R)^2 == ", format(round(r2, 2), nsmall = 2))
N_diet <- ggplot(data_summary, aes(x = FoodTroph, y = n)) +
  geom_segment(aes(xend = FoodTroph, yend = Fitted), linetype = "dotted") +
  geom_point(aes(color = Above, fill = Color), shape = 21, size = 3, stroke = 1, show.legend = F) +
  scale_fill_identity() +
  scale_color_manual(values = c("above" = "black", "below" = "black")) +
  stat_smooth(method = "lm", formula = y ~ I(x - 1) + 0, color = "mediumpurple", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = r2_label, parse = TRUE,
           hjust = 1.1, vjust = 1.5, size = 5) + theme_classic() +
  scale_x_continuous(name = "Trophic level", limits = c(2, 4.5)) +
  scale_y_continuous(name = "N diet (%)", limits = c(2, 12), breaks = seq(2,10,2)) +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

model_P      <- lm(p ~ I(FoodTroph - 1) + 0, data = data_summary)
data_summary <- data_summary |> mutate(Fitted = predict(model_P),
    Above = ifelse(p > Fitted, "above", "below"), Color = ifelse(Above == "above", "dodgerblue2", "coral2"))
r2           <- summary(model_P)$r.squared
r2_label     <- paste0("italic(R)^2 == ", format(round(r2, 2), nsmall = 2))
P_diet <- ggplot(data_summary, aes(x = FoodTroph, y = p)) +
  geom_segment(aes(xend = FoodTroph, yend = Fitted), linetype = "dotted") +
  geom_point(aes(color = Above, fill = Color), shape = 21, size = 3, stroke = 1, show.legend = F) +
  scale_fill_identity() +
  scale_color_manual(values = c("above" = "black", "below" = "black")) +
  stat_smooth(method = "lm", formula = y ~ I(x - 1) + 0, color = "mediumpurple", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = r2_label, parse = TRUE,
           hjust = 1.1, vjust = 1.5, size = 5) + theme_classic() +
  scale_x_continuous(name = "Trophic level", limits = c(2, 4.5)) +
  scale_y_continuous(name = "P diet (%)") +
  theme(axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

Figure_S2 = C_diet + N_diet + P_diet + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Add estimates
nutrients_imputed = nutrients_imputed |>
  mutate(Species = gsub("_", " ", Species)) |>
  left_join((rfishbase::load_taxa() |> select(SpecCode, Species)), by = "Species") |> 
  left_join((rfishbase::ecology() |> select(SpecCode, FoodTroph) |> group_by(SpecCode) |> 
               summarise(FoodTroph = mean(FoodTroph))), by = "SpecCode") |> 
  dplyr::select(-SpecCode) |> 
  mutate(Diet_C = predict(model_C, newdata = pick(everything())),
         Diet_N = predict(model_N, newdata = pick(everything())),
         Diet_P = predict(model_P, newdata = pick(everything())))

#### 3. Constantes       ----

aC   = 0.8
aN   = 0.8
aP   = 0.7
F0Nz = 3.7e-03
F0Pz = 3.7e-04

#### 4. Growth           ----
lw_growth <- rfishbase::length_weight() |> group_by(SpecCode) |>
  summarise(lwa = mean(a, na.rm = T), lwb = mean(b, na.rm = T)) |>
  left_join(rfishbase::popgrowth() |> group_by(SpecCode) |>
    summarise(linf = mean(Loo, na.rm = T), K = mean(K, na.rm = T), t0 = mean(to, na.rm = T)), by = "SpecCode") |>
  left_join(rfishbase::load_taxa() |> select(SpecCode, Species),  by = "SpecCode") |> select(Species, everything()) |>
  mutate(t0 = ifelse(is.nan(t0), NA, t0)) |> dplyr::select(-SpecCode)

lw_growth$Species <- gsub(" ", "_", lw_growth$Species)
phy_lw            <- fishtree_phylogeny(species = lw_growth$Species)
phy_lw$tip.label  <- gsub(" ", "_", phy_lw$tip.label)
lw_growth         <- lw_growth |> filter(!is.na(Species))

lw_growth         <- lw_growth |>
  impute_trait_phylopars(trait_col = "lwa", phy = phy_lw) |>
  impute_trait_phylopars(trait_col = "lwb", phy = phy_lw) |>
  impute_trait_phylopars(trait_col = "linf", phy = phy_lw) |> 
  impute_trait_phylopars(trait_col = "K", phy = phy_lw) |>
  impute_trait_phylopars(trait_col = "t0", phy = phy_lw) 

#### 5. Caudal fin       ----
Caudal_fin_FB <- rfishbase::morphometrics() |> data.frame() |> 
  dplyr::select(SpecCode, AspectRatio) |> 
  left_join(rfishbase::load_taxa(), by = "SpecCode") |> 
  dplyr::select(Species, Family, AspectRatio) |> 
  group_by(Species, Family) |> summarise(AspectRatio = mean(AspectRatio, na.rm = TRUE)) |>
  mutate(AspectRatio = ifelse(is.nan(AspectRatio), NA, AspectRatio)) |> ungroup()

Caudal_fin_FB$Species <- gsub(" ", "_", Caudal_fin_FB$Species)
phy_caudal <- fishtree_phylogeny(species = Caudal_fin_FB$Species)
phy_caudal$tip.label <- gsub(" ", "_", phy_caudal$tip.label)
caudal_imputed <- Caudal_fin_FB |> filter(!is.na(Species)) |> 
  impute_trait_phylopars(trait_col = "AspectRatio", phy = phy_caudal)

#### 6. Metabolism       ----
## Looking for f0 and alpha numerically
# Dataset from Rosen et al. 2025
Metabolic_Rosen = Metabolic_Rosen |> 
  mutate(SMR_gC_day = SMR * 0.312 * 24 / 1000, MMR_gC_day = MMR * 0.312 * 24 / 1000,
         log_SMR_gC_day = log(SMR_gC_day), log_MMR_gC_day = log(MMR_gC_day),
         log_weight = log(Weight),
         invKT = 1 / (8.617e-5 * (Temperature + 273.15))) 
# Dataset from Schiettekatte et al. 2021
Metabolic_Nina = Metabolic_Nina |> 
  rename(SMR_gC_day = SMR, MMR_gC_day = MaxMR) |> 
  mutate(SMR_gC_day = SMR_gC_day * 0.312 * 24 / 1000, MMR_gC_day = MMR_gC_day * 0.312 * 24 / 1000,
         log_SMR_gC_day = log(SMR_gC_day), log_MMR_gC_day = log(MMR_gC_day),
         log_weight = log(Weight..kg. * 1000),
         invKT = 1 / (8.617e-5 * (MeanTemp...C. + 273.15))) |> 
  dplyr::filter(Etat == "ok", SMR_gC_day > 0)

# Combine datasets
Metabolic = data.frame(SMR_gC_day = c(Metabolic_Rosen$SMR_gC_day, Metabolic_Nina$SMR_gC_day),
                       MMR_gC_day = c(Metabolic_Rosen$MMR_gC_day, Metabolic_Nina$MMR_gC_day),
                       log_SMR_gC_day = log(c(Metabolic_Rosen$SMR_gC_day, Metabolic_Nina$SMR_gC_day)),
                       log_MMR_gC_day = log(c(Metabolic_Rosen$MMR_gC_day, Metabolic_Nina$MMR_gC_day)),
                       log_weight = log(c(Metabolic_Rosen$Weight, Metabolic_Nina$Weight..kg.*1000)),
                       Weight = c(Metabolic_Rosen$Weight, Metabolic_Nina$Weight..kg.*1000),
                       Species = c(Metabolic_Rosen$Species, Metabolic_Nina$Species),
                       invKT = c(Metabolic_Rosen$invKT, Metabolic_Nina$invKT),
                       Temperature = c(Metabolic_Rosen$Temperature, Metabolic_Nina$MeanTemp...C.),
                       Dataset = c(rep("Rosen", length(Metabolic_Rosen$invKT)), 
                                   rep("Nina", length(Metabolic_Nina$invKT))))

# Test to explore mathematically alpha, E, f0 and theta
(alpha <- coef(lm(log_SMR_gC_day ~ log_weight + invKT, data = Metabolic))["log_weight"])
(E     <- coef(lm(log_SMR_gC_day ~ log_weight + invKT, data = Metabolic))["invKT"])
Temp_ref <- 27 # Mainly tropical reef fishes
Metabolic_summary = Metabolic |> 
  mutate(f0 = SMR_gC_day / (Weight^alpha * exp(E * (1 / (Temp_ref + 273.15) - 1 / (Temperature + 273.15)) / 8.617e-5)),
         theta = (SMR_gC_day + MMR_gC_day) / (2 * SMR_gC_day)) |> 
  group_by(Species) |> 
  summarise(f0_avg = mean(f0), theta_avg = mean(theta), f0_sd = sd(f0), theta_sd = sd(theta))

# Perform relationship for SMR
SMR_Weight_relationship <- brm(log_SMR_gC_day ~ log_weight, data = Metabolic,
  family = gaussian(), prior = c(prior(normal(0, 5), class = "Intercept"), prior(normal(0.75, 0.5), class = "b")),
  chains = 4, cores = 4, iter = 4000, seed = 123)
Metabolic_model_data          <- model.frame(SMR_Weight_relationship)
Metabolic_model_data$SMR_pred <- exp(fitted(SMR_Weight_relationship, scale = "linear")[, "Estimate"])
Metabolic_model_data$SMR_obs  <- exp(Metabolic_model_data$log_SMR_gC_day)
Metabolic_model_data$Weight   <- exp(Metabolic_model_data$log_weight)

slope     <- round(fixef(SMR_Weight_relationship)["log_weight", "Estimate"], 3)
intercept <- round(exp(fixef(SMR_Weight_relationship)["Intercept", "Estimate"]), 3)
r2_bayes  <- round(bayes_R2(SMR_Weight_relationship)[1], 2)
eq_label  <- bquote(italic(SMR) == .(intercept) %*% Weight^.(slope))
r2_label  <- bquote(italic(R)^2 == .(r2_bayes))

label_data_SMR <- data.frame(
  x = 0.005,
  y_eq = 1 * 0.8,
  y_r2 = 1 * 0.7, 
  eq_label = as.character(deparse(eq_label)),
  r2_label = as.character(deparse(r2_label)))

Figure_S3_A <- Metabolic_model_data |>
  mutate(Above = ifelse(SMR_obs > SMR_pred, "above", "below"),
         Color = ifelse(Above == "above", "dodgerblue2", "coral2")) |> 
  ggplot(aes(x = Weight, y = SMR_obs)) +
  geom_segment(aes(xend = Weight, yend = SMR_pred), linetype = "dotted") +
  geom_point(aes(fill = Color, color = Above), shape = 21, size = 3, stroke = 1, show.legend = FALSE) +
  scale_fill_identity() + theme_classic() +
  scale_color_manual(values = c("above" = "black", "below" = "black")) +
  geom_line(aes(y = SMR_pred), color = "mediumpurple", linewidth = 1.2) +
  scale_x_log10(name = "Weight (g)", labels = label_number(accuracy = 0.01),
                limits = c(0.001, 1000)) +
  scale_y_log10(name = expression("Standard Metabolic Rate"~(g~C~day^{-1})), labels = label_number(accuracy = 0.0001),
                limits = c(0.0001, 1)) +
  geom_text(data = label_data_SMR, aes(x = x, y = y_eq, label = eq_label), 
            hjust = 0, vjust = 0, size = 5, parse = TRUE) +
  geom_text(data = label_data_SMR, aes(x = x, y = y_r2, label = r2_label), 
            hjust = 0, vjust = 1.5, size = 5, parse = TRUE) +
  theme(axis.title = element_text(size = 16), axis.text  = element_text(size = 14))

# Perform relationship for MMR
MMR_Weight_relationship <- brm(log_MMR_gC_day ~ log_weight, data = Metabolic,
  family = gaussian(), prior = c(prior(normal(0, 5), class = "Intercept"), prior(normal(0.75, 0.5), class = "b")),
  chains = 4, cores = 4, iter = 4000, seed = 123)
Metabolic_model_data          <- model.frame(MMR_Weight_relationship)
Metabolic_model_data$MMR_pred <- exp(fitted(MMR_Weight_relationship, scale = "linear")[, "Estimate"])
Metabolic_model_data$MMR_obs  <- exp(Metabolic_model_data$log_MMR_gC_day)
Metabolic_model_data$Weight   <- exp(Metabolic_model_data$log_weight)

slope     <- round(fixef(MMR_Weight_relationship)["log_weight", "Estimate"], 3)
intercept <- round(exp(fixef(MMR_Weight_relationship)["Intercept", "Estimate"]), 3)
r2_bayes  <- round(bayes_R2(MMR_Weight_relationship)[1], 2)
eq_label  <- bquote(italic(MMR) == .(intercept) %*% Weight^.(slope))
r2_label  <- bquote(italic(R)^2 == .(r2_bayes))

label_data_MMR <- data.frame(
  x = 0.005,
  y_eq = 1 * 0.8,
  y_r2 = 1 * 0.7, 
  eq_label = as.character(deparse(eq_label)),
  r2_label = as.character(deparse(r2_label)))

Figure_S3_B <- Metabolic_model_data |>
  mutate(Above = ifelse(MMR_obs > MMR_pred, "above", "below"),
         Color = ifelse(Above == "above", "dodgerblue2", "coral2")) |> 
  ggplot(aes(x = Weight, y = MMR_obs)) +
  geom_segment(aes(xend = Weight, yend = MMR_pred), linetype = "dotted") +
  geom_point(aes(fill = Color, color = Above), shape = 21, size = 3, stroke = 1, show.legend = FALSE) +
  scale_fill_identity() + theme_classic() +
  scale_color_manual(values = c("above" = "black", "below" = "black")) +
  geom_line(aes(y = MMR_pred), color = "mediumpurple", linewidth = 1.2) +
  scale_x_log10(name = "Weight (g)", labels = label_number(accuracy = 0.01),
                limits = c(0.001, 1000)) +
  scale_y_log10(name = expression("Maximum Metabolic Rate"~(g~C~day^{-1})), labels = label_number(accuracy = 0.0001),
                limits = c(0.0001, 1)) +
  geom_text(data = label_data_MMR, aes(x = x, y = y_eq, label = eq_label), 
            hjust = 0, vjust = 0, size = 5, parse = TRUE) +
  geom_text(data = label_data_MMR, aes(x = x, y = y_r2, label = r2_label), 
            hjust = 0, vjust = 1.5, size = 5, parse = TRUE) +
  theme(axis.title = element_text(size = 16), axis.text  = element_text(size = 14))

(Figure_S3 = Figure_S3_A + Figure_S3_B + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16), axis.text  = element_text(size = 14)))

#### 7. Compilation  ----

### Clean Species info
# Fix spp. rows
Western_Med <- B_Useful_species_list |>
  distinct(Species, Genus, `C/DM`, `N/DM`, `P/DM`) |>
  mutate(Genus = ifelse(grepl("spp\\.", Species), word(Species, 1), Genus)) |>
  left_join(load_taxa() |> select(Genus, Family) |> distinct(), by = "Genus") |>
  relocate(Family, .after = Genus)
# Fix names
Western_Med$Species[Western_Med$Species == "Trigloporus lastoviza"] = "Chelidonichthys lastoviza"
Western_Med$Species[Western_Med$Species == "Pteromylaeus bovinus"]  = "Aetomylaeus bovinus"
Western_Med$Species[Western_Med$Species == "Liza ramada"]           = "Chelon ramada"
Western_Med$Species[Western_Med$Species == "Liza aurata"]           = "Chelon auratus"
Western_Med$Species[Western_Med$Species == "Liza saliens"]          = "Chelon saliens"
# Fix non-found species
Western_Med = Western_Med |> mutate(Genus = ifelse(is.na(Genus) & is.na(Family) | Genus == "", word(Species, 1), Genus),
                             Genus = ifelse(grepl("spp\\.", Species), word(Species, 1), Genus)) |>
  left_join(load_taxa() |> select(Genus, Family) |> distinct(), by = "Genus") |> 
  relocate(Family.y, .after = Genus) |> dplyr::select(-Family.x) |> rename(Family = Family.y)
# Fix Family rows
Western_Med$Family[Western_Med$Species == "Myctophidae"] = "Myctophidae"  
Western_Med$Genus[Western_Med$Species  == "Myctophidae"] = NA
Western_Med$Family[Western_Med$Species == "Blenniidae"]  = "Blenniidae"
Western_Med$Genus[Western_Med$Species  == "Blenniidae"]  = NA
# Formating the data
Western_Med = Western_Med |> rename(Qc = `C/DM`, Qn = `N/DM`, Qp = `P/DM`) |> 
  mutate(Dataset = "B_Useful_MED_W")

### Add CNP body composition
Nut_West = nutrients_imputed |> select(Species, Family, CDM, NDM, PDM) |> 
  rename(Qc = CDM, Qn = NDM, Qp = PDM) |> mutate(Dataset = "Schiettekatte_2021") |> 
  left_join(load_taxa() |> select(Species, Genus) |> distinct(), by = "Species") |> 
  relocate(Genus, .before = Family)
Western_Med = bind_rows(Western_Med, Nut_West)
# Phylogenetic workflow
Western_Med$Species <- gsub(" ", "_", Western_Med$Species)
phy_lw              <- fishtree_phylogeny(species = Western_Med$Species)
phy_lw$tip.label    <- gsub(" ", "_", phy_lw$tip.label)
# Body composition
Western_Med         <- Western_Med |> filter(!is.na(Species)) |> 
  impute_trait_phylopars(trait_col = "Qc", phy = phy_lw) |>
  impute_trait_phylopars(trait_col = "Qn", phy = phy_lw) |>
  impute_trait_phylopars(trait_col = "Qp", phy = phy_lw) |> 
  select(-c(Qc_type, Qn_type, Qp_type)) |> relocate(Dataset, .after = Qp) |> group_by(Species) |>
  slice_max(order_by = (Dataset == "B_Useful_MED_W"), with_ties = FALSE) |> ungroup() |> 
  fill_trait_hierarchy("Qc") %>%
  fill_trait_hierarchy("Qn") %>%
  fill_trait_hierarchy("Qp")  

### Add CNP diet
Western_Med         <- Western_Med |> mutate(Species = gsub("_", " ", Species)) |> 
  left_join(rfishbase::load_taxa(), by = "Species") |> 
  left_join(rfishbase::ecology(), by = "SpecCode", relationship = "many-to-many") |> 
  select(c(1:7,51)) |> mutate(Species = gsub(" ", "_", Species)) |> 
  rename(Genus = Genus.x, Family = Family.x) |> 
  impute_trait_phylopars(trait_col = "FoodTroph", phy = phy_lw) |> 
  fill_trait_hierarchy("FoodTroph") |> select(-FoodTroph_type) |> 
  mutate(Dc = predict(model_C, newdata = pick(everything())),
         Dn = predict(model_N, newdata = pick(everything())),
         Dp = predict(model_P, newdata = pick(everything()))) |> 
  relocate(c(FoodTroph, Dataset), .after = Dp) |> mutate(Qc = Qc * 100, Qn = Qn * 100, Qp = Qp * 100) |> 
  dplyr::filter(Dataset == "B_Useful_MED_W") |> rename(h = FoodTroph)

### Add constantes
Western_Med         <- Western_Med |> mutate(ac = 0.8, an = 0.8, ap = 0.7, f0nz = 3.7e-03, f0pz = 3.7e-04) |> 
  relocate(Dataset, .after = f0pz)

#### 8. Export the data  ----
ggsave(Figure_S1, filename = "Figure_S1.png", path = "Outputs/", device = "png", width = 4,  height = 7.5, dpi = 300)  
ggsave(Figure_S2, filename = "Figure_S2.png", path = "Outputs/", device = "png", width = 10, height = 3.5, dpi = 300) 
ggsave(Figure_S3, filename = "Figure_S3.png", path = "Outputs/", device = "png", width = 10, height = 5.0, dpi = 300) 
