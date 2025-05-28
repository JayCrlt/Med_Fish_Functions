rm(list = ls())
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl')
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars")
library("ggridges") ; library("patchwork") ; library("fishflux")

## What is a function? A function is a rate. (Jax)
## Then, the idea would be to define three core functions: i) Biomass production through Growth rates;
## ii) Predation rate through  Q/B * Biomass; iii) nutrient cycling through C,N,P/DM + FishFlux
## All those data are available in FishBase, but maybe not directly from RFishBase, then, we would have to compile 
## The data by ourselve

# Raw data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Medits_total       <- read.delim("Data/TATB_WMED_1999-2021_clean.csv", sep = ";")
Nutrients_FishBase <- readxl::read_excel("Data/Nutrients_FishBase.xlsx") #87 values for C, 206 for N, 46 for P
# Nina data CNP diet
cnp_diet           <- read.delim("Data/cnp_diet.csv", sep = ",") |> 
  rename(Species = species) |> 
  left_join(rfishbase::load_taxa(), by = "Species") |> 
  left_join(rfishbase::ecology(), by = "SpecCode") 

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

#### 1. CNP Body mass
# Distinct species in Nutrients database
Nutrients_FishBase = Nutrients_FishBase |> group_by(Species, Family) |> 
  summarise(`C/DM` = mean(`C/DM`), `N/DM` = mean(`N/DM`), `P/DM` = mean(`P/DM`))
# Standardize species names once
Nutrients_FishBase$Species <- gsub(" ", "_", Nutrients_FishBase$Species)
# Load phylogeny once
phy <- fishtree_phylogeny(species = Nutrients_FishBase$Species)
phy$tip.label <- gsub(" ", "_", phy$tip.label)
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

#### 2. CNP diets
data_summary <- cnp_diet |> group_by(Species) |> summarise_all(mean, na.rm = TRUE) |> filter(!is.na(FoodTroph), !is.na(c))
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

#### 0. Export the data
