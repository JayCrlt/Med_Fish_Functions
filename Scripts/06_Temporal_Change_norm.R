#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom") ; library("gstat") ; library("forcats")

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

### Define the Fish hexagonal grid
Hexagonal_grid_fish_surveys <- Hexagonal_grid |> 
  dplyr::filter(id %in% (Medit_FunCatch_without_NA |> dplyr::distinct(HEX_ID) |> dplyr::pull(HEX_ID)) | 
                  grid_id %in% (Medit_FunCatch_without_NA |> dplyr::distinct(HEX_ID) |> dplyr::pull(HEX_ID)))

### Define the fish dataset for the last five years
Medit_FunCatch_last_five_years = Medit_FunCatch_without_NA |> dplyr::filter(YEAR %in% c(2019:2021))

### Define the Fish hexagonal grid for the last five years to compare
Hexagonal_grid_five_years <- Hexagonal_grid |> 
  dplyr::filter(id %in% (Medit_FunCatch_last_five_years |> dplyr::distinct(HEX_ID) |> dplyr::pull(HEX_ID)) | 
                  grid_id %in% (Medit_FunCatch_last_five_years |> dplyr::distinct(HEX_ID) |> dplyr::pull(HEX_ID)))

### Keep the highest values for the 5 functions + multifunctional index from the last three years
medits_sf_percentile_hex <- st_join(medits_sf_percentile, Hexagonal_grid |> dplyr::select(id, grid_id)) |> 
  mutate(HEX_ID = coalesce(id, grid_id)) |> dplyr::select(-c(id, grid_id)) |> dplyr::filter(YEAR %in% c(2019:2021)) |> 
  st_drop_geometry() |> group_by(HEX_ID) |> summarise(max_Fn  = max(community_Fn, na.rm = TRUE), 
                                                        max_Fp  = max(community_Fp, na.rm = TRUE),
                                                        max_Gc  = max(community_Gc, na.rm = TRUE),
                                                        max_IcP = max(Ic_plank, na.rm = TRUE),
                                                        max_IcB = max(Ic_benthivorous, na.rm = TRUE),
                                                        max_MF  = max(Multifunctionality, na.rm = TRUE)) |> ungroup() |>
  left_join(Hexagonal_grid |> select(id, grid_id, geometry), by = c("HEX_ID" = "id")) |> st_as_sf()

### See if the spatial interpolation is reasonable
ggplot() +
  geom_sf(data = Hexagonal_grid_fish_surveys, fill = "skyblue", color = "blue", size = 0.2, alpha = 0.5) + 
  geom_sf(data = medits_sf_percentile_hex, fill = "red", color = "orange", size = 0.3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Comparison of Hexagonal Grids", subtitle = "Fish surveys all years (blue), Last 3 years (red)") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  theme(panel.grid = element_blank())

### Spatial interpolation
interp_points <- st_centroid(medits_sf_percentile_hex)
interp_grid   <- st_centroid(Hexagonal_grid_fish_surveys)
interp_Fn  <- predict(gstat(formula = max_Fn  ~ 1, locations = interp_points), newdata = interp_grid)
interp_Fp  <- predict(gstat(formula = max_Fp  ~ 1, locations = interp_points), newdata = interp_grid)
interp_Gc  <- predict(gstat(formula = max_Gc  ~ 1, locations = interp_points), newdata = interp_grid)
interp_IcP <- predict(gstat(formula = max_IcP ~ 1, locations = interp_points), newdata = interp_grid)
interp_IcB <- predict(gstat(formula = max_IcB ~ 1, locations = interp_points), newdata = interp_grid)
interp_MF  <- predict(gstat(formula = max_MF  ~ 1, locations = interp_points), newdata = interp_grid)
Hexagonal_grid_fish_surveys_max = Hexagonal_grid_fish_surveys |> 
  mutate(interp_Fn_max = interp_Fn$var1.pred, interp_Fp_max = interp_Fp$var1.pred, interp_Gc_max = interp_Gc$var1.pred,
         interp_IcP_max = interp_IcP$var1.pred, interp_IcB_max = interp_IcB$var1.pred, 
         interp_Mf_max = interp_MF$var1.pred) |> 
  rename(HEX_ID = id)

### Keep the highest values for the 5 functions + multifunctional index prior to 2019
medits_sf_percentile_hex_prev <- st_join(medits_sf_percentile, Hexagonal_grid |> dplyr::select(id, grid_id)) |> 
  mutate(HEX_ID = coalesce(id, grid_id)) |> dplyr::select(-c(id, grid_id)) |> dplyr::filter(YEAR %notin% c(2019:2021)) |> 
  st_drop_geometry() |> group_by(HEX_ID) |> summarise(max_Fn  = max(community_Fn, na.rm = TRUE), 
                                                        max_Fp  = max(community_Fp, na.rm = TRUE),
                                                        max_Gc  = max(community_Gc, na.rm = TRUE),
                                                        max_IcP = max(Ic_plank, na.rm = TRUE),
                                                        max_IcB = max(Ic_benthivorous, na.rm = TRUE),
                                                        max_MF  = max(Multifunctionality, na.rm = TRUE)) |> ungroup() |>
  left_join(Hexagonal_grid |> select(id, grid_id, geometry), by = c("HEX_ID" = "id")) |> st_as_sf()

### Compare to previous years
comp_hex <- Hexagonal_grid_fish_surveys_max |>
  left_join(medits_sf_percentile_hex_prev |> st_drop_geometry() |>  
              select(HEX_ID, max_Fn, max_Fp, max_Gc, max_IcP, max_IcB, max_MF), by = "HEX_ID") |>
  mutate(diff_Fn= interp_Fn_max - max_Fn, diff_Fp = interp_Fp_max - max_Fp, diff_Gc = interp_Gc_max - max_Gc,
    diff_IcP = interp_IcP_max - max_IcP, diff_IcB = interp_IcB_max - max_IcB, diff_MF = interp_Mf_max - max_MF) |>
  mutate(diff_Fn_cat  = case_when(diff_Fn  > 0 ~ "Positive", diff_Fn  < 0 ~ "Negative", TRUE ~ "NA"),
         diff_Fp_cat  = case_when(diff_Fp  > 0 ~ "Positive", diff_Fp  < 0 ~ "Negative", TRUE ~ "NA"),
         diff_Gc_cat  = case_when(diff_Gc  > 0 ~ "Positive", diff_Gc  < 0 ~ "Negative", TRUE ~ "NA"),
         diff_IcP_cat = case_when(diff_IcP > 0 ~ "Positive", diff_IcP < 0 ~ "Negative", TRUE ~ "NA"),
         diff_IcB_cat = case_when(diff_IcB > 0 ~ "Positive", diff_IcB < 0 ~ "Negative", TRUE ~ "NA"),
         diff_MF_cat  = case_when(diff_MF  > 0 ~ "Positive", diff_MF  < 0 ~ "Negative", TRUE ~ "NA")) |> drop_na() |>
  st_drop_geometry() |> st_as_sf(coords = c("grid_lon", "grid_lat"), crs = 4326)

### Plots
label_val_Gc <- table(comp_hex$diff_Gc_cat) |> as.data.frame() |> 
  mutate(Freq = paste0(round(Freq / 2487 * 100, 2), "%")) |> dplyr::filter(Var1 == "Positive") |> pull(Freq)
Figure_S3A <- ggplot() +
  geom_sf(data = subset(comp_hex, diff_Gc_cat == "Negative"),
          aes(fill = diff_Gc_cat, color = diff_Gc_cat, size = diff_Gc_cat, shape = diff_Gc_cat)) +
  geom_sf(data = subset(comp_hex, diff_Gc_cat == "Positive"),
          aes(fill = diff_Gc_cat, color = diff_Gc_cat, size = diff_Gc_cat, shape = diff_Gc_cat)) +
  scale_shape_manual(values = c("Positive" = 21, "Negative" = 20)) +
  scale_fill_manual( values = c("Positive" = "#B4CBF0", "Negative" = "grey")) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "grey50")) +
  scale_size_manual( values = c("Positive" = 2, "Negative" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_label(aes(x = 32, y = 45, label = label_val_Gc), inherit.aes = F, size = 5, fill = "#B4CBF0", color = "white") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Production")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 0),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

label_val_Fn <- table(comp_hex$diff_Fn_cat) |> as.data.frame() |> 
  mutate(Freq = paste0(round(Freq / 2487 * 100, 2), "%")) |> dplyr::filter(Var1 == "Positive") |> pull(Freq)
Figure_S3B <- ggplot() +
  geom_sf(data = subset(comp_hex, diff_Fn_cat == "Negative"),
          aes(fill = diff_Fn_cat, color = diff_Fn_cat, size = diff_Fn_cat, shape = diff_Fn_cat)) +
  geom_sf(data = subset(comp_hex, diff_Fn_cat == "Positive"),
          aes(fill = diff_Fn_cat, color = diff_Fn_cat, size = diff_Fn_cat, shape = diff_Fn_cat)) +
  scale_shape_manual(values = c("Positive" = 21, "Negative" = 20)) +
  scale_fill_manual( values = c("Positive" = "#A3B79C", "Negative" = "grey")) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "grey50")) +
  scale_size_manual( values = c("Positive" = 2, "Negative" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_label(aes(x = 32, y = 45, label = label_val_Fn), inherit.aes = F, size = 5, fill = "#A3B79C", color = "white") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Nitrogen excretion")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 0),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

label_val_Fp <- table(comp_hex$diff_Fp_cat) |> as.data.frame() |> 
  mutate(Freq = paste0(round(Freq / 2487 * 100, 2), "%")) |> dplyr::filter(Var1 == "Positive") |> pull(Freq)
Figure_S3C <- ggplot() +
  geom_sf(data = subset(comp_hex, diff_Fp_cat == "Negative"),
          aes(fill = diff_Fp_cat, color = diff_Fp_cat, size = diff_Fp_cat, shape = diff_Fp_cat)) +
  geom_sf(data = subset(comp_hex, diff_Fp_cat == "Positive"),
          aes(fill = diff_Fp_cat, color = diff_Fp_cat, size = diff_Fp_cat, shape = diff_Fp_cat)) +
  scale_shape_manual(values = c("Positive" = 21, "Negative" = 20)) +
  scale_fill_manual( values = c("Positive" = "#F2DE93", "Negative" = "grey")) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "grey50")) +
  scale_size_manual( values = c("Positive" = 2, "Negative" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_label(aes(x = 32, y = 45, label = "47.00%"), inherit.aes = F, size = 5, fill = "#F2DE93", color = "white") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Phosphorus excretion")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 0),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

label_val_IcP <- table(comp_hex$diff_IcP_cat) |> as.data.frame() |> 
  mutate(Freq = paste0(round(Freq / 2487 * 100, 2), "%")) |> dplyr::filter(Var1 == "Positive") |> pull(Freq)
Figure_S3D <- ggplot() +
  geom_sf(data = subset(comp_hex, diff_IcP_cat == "Negative"),
          aes(fill = diff_IcP_cat, color = diff_IcP_cat, size = diff_IcP_cat, shape = diff_IcP_cat)) +
  geom_sf(data = subset(comp_hex, diff_IcP_cat == "Positive"),
          aes(fill = diff_IcP_cat, color = diff_IcP_cat, size = diff_IcP_cat, shape = diff_IcP_cat)) +
  scale_shape_manual(values = c("Positive" = 21, "Negative" = 20)) +
  scale_fill_manual( values = c("Positive" = "#CCA9DD", "Negative" = "grey")) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "grey50")) +
  scale_size_manual( values = c("Positive" = 2, "Negative" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_label(aes(x = 32, y = 45, label = label_val_IcP), inherit.aes = F, size = 5, fill = "#CCA9DD", color = "white") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Planktivory")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 0),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

label_val_IcB <- table(comp_hex$diff_IcB_cat) |> as.data.frame() |> 
  mutate(Freq = paste0(round(Freq / 2487 * 100, 2), "%")) |> dplyr::filter(Var1 == "Positive") |> pull(Freq)
Figure_S3E <- ggplot() +
  geom_sf(data = subset(comp_hex, diff_IcB_cat == "Negative"),
          aes(fill = diff_IcB_cat, color = diff_IcB_cat, size = diff_IcB_cat, shape = diff_IcB_cat)) +
  geom_sf(data = subset(comp_hex, diff_IcB_cat == "Positive"),
          aes(fill = diff_IcB_cat, color = diff_IcB_cat, size = diff_IcB_cat, shape = diff_IcB_cat)) +
  scale_shape_manual(values = c("Positive" = 21, "Negative" = 20)) +
  scale_fill_manual( values = c("Positive" = "#FAC898", "Negative" = "grey")) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "grey50")) +
  scale_size_manual( values = c("Positive" = 2, "Negative" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_label(aes(x = 32, y = 45, label = label_val_IcB), inherit.aes = F, size = 5, fill = "#FAC898", color = "white") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Benthivory")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 0),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

label_val_Mf <- table(comp_hex$diff_MF_cat) |> as.data.frame() |> 
  mutate(Freq = paste0(round(Freq / 2487 * 100, 2), "%")) |> dplyr::filter(Var1 == "Positive") |> pull(Freq)
Figure_S3F <- ggplot() +
  geom_sf(data = subset(comp_hex, diff_MF_cat == "Negative"),
          aes(fill = diff_MF_cat, color = diff_MF_cat, size = diff_MF_cat, shape = diff_MF_cat)) +
  geom_sf(data = subset(comp_hex, diff_MF_cat == "Positive"),
          aes(fill = diff_MF_cat, color = diff_MF_cat, size = diff_MF_cat, shape = diff_MF_cat)) +
  scale_shape_manual(values = c("Positive" = 21, "Negative" = 20)) +
  scale_fill_manual( values = c("Positive" = "#FF968A", "Negative" = "grey")) +
  scale_color_manual(values = c("Positive" = "black", "Negative" = "grey50")) +
  scale_size_manual( values = c("Positive" = 2, "Negative" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  geom_label(aes(x = 32, y = 45, label = label_val_Mf), inherit.aes = F, size = 5, fill = "#FF968A", color = "white") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Multifunctional Index")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 0),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

# Altogether
Figure_S3 = Figure_S3A + Figure_S3B + Figure_S3C + Figure_S3D + Figure_S3E + Figure_S3F + 
  plot_layout(guides = "collect", ncol = 2) 

### Ridgelplots of functions over time
ridge_summary <- medits_sf_percentile %>%
  st_drop_geometry() %>% mutate(period = case_when(YEAR %in% 1999:2001 ~ "1999-2001", YEAR %in% 2002:2006 ~ "2002-2006",
    YEAR %in% 2007:2011 ~ "2007-2011", YEAR %in% 2012:2016 ~ "2012-2016", YEAR %in% 2017:2021 ~ "2017-2021",
    TRUE ~ NA_character_)) %>%filter(!is.na(period)) %>%
  select(period, community_Fn, community_Fp, community_Gc, Ic_plank, Ic_benthivorous, Multifunctionality) %>%
  pivot_longer(cols = -period, names_to = "Function", values_to = "Value") %>% group_by(period, Function) %>%
  summarise(mean = median(Value, na.rm = TRUE), q25  = quantile(Value, 0.25, na.rm = TRUE), 
            q75  = quantile(Value, 0.75, na.rm = TRUE), q05  = quantile(Value, 0.05, na.rm = TRUE),
            q95  = quantile(Value, 0.95, na.rm = TRUE), .groups = "drop")

Figure_S4A <- ridge_summary %>% dplyr::filter(Function == "community_Gc") |> 
  ggplot(aes(y = fct_rev(factor(period)), x = mean)) +
  geom_linerange(aes(xmin = q05, xmax = q95), color = "#cfe2f3", size = 1) +   scale_x_log10() +
  geom_linerange(aes(xmin = q25, xmax = q75), color = "#B4CBF0", size = 2) + 
  geom_point(fill = "#B4CBF0", size = 5, shape = 21, color = "black") +
  labs(y = "", x = expression("Production ("*gC~m^{-2}~d^{-1}*")")) + theme_minimal(base_size = 14) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_S4B <- ridge_summary %>% dplyr::filter(Function == "community_Fn") |> 
  ggplot(aes(y = fct_rev(factor(period)), x = mean)) +
  geom_linerange(aes(xmin = q05, xmax = q95), color = "#d9ead3", size = 1) +   scale_x_log10() +
  geom_linerange(aes(xmin = q25, xmax = q75), color = "#A3B79C", size = 2) + 
  geom_point(fill = "#A3B79C", size = 5, shape = 21, color = "black") +
  labs(y = "", x = expression("Nitrogen excretion ("*gN~m^{-2}~d^{-1}*")")) + theme_minimal(base_size = 14) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_S4C <- ridge_summary %>% dplyr::filter(Function == "community_Fp") |> 
  ggplot(aes(y = fct_rev(factor(period)), x = mean)) +
  geom_linerange(aes(xmin = q05, xmax = q95), color = "#fff2cc", size = 1) +   scale_x_log10() +
  geom_linerange(aes(xmin = q25, xmax = q75), color = "#F2DE93", size = 2) + 
  geom_point(fill = "#F2DE93", size = 5, shape = 21, color = "black") +
  labs(y = "", x = expression("Phosphorus excretion ("*gP~m^{-2}~d^{-1}*")")) + theme_minimal(base_size = 14) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_S4D <- ridge_summary %>% dplyr::filter(Function == "Ic_plank") |> 
  ggplot(aes(y = fct_rev(factor(period)), x = mean)) +
  geom_linerange(aes(xmin = q05, xmax = q95), color = "#d9d2e9", size = 1) +   scale_x_log10() +
  geom_linerange(aes(xmin = q25, xmax = q75), color = "#CCA9DD", size = 2) + 
  geom_point(fill = "#CCA9DD", size = 5, shape = 21, color = "black") +
  labs(y = "", x = expression("Planktivory ("*gC~m^{-2}~d^{-1}*")")) + theme_minimal(base_size = 14) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_S4E <- ridge_summary %>% dplyr::filter(Function == "Ic_benthivorous") |> 
  ggplot(aes(y = fct_rev(factor(period)), x = mean)) +
  geom_linerange(aes(xmin = q05, xmax = q95), color = "#fff2cc", size = 1) +   scale_x_log10() +
  geom_linerange(aes(xmin = q25, xmax = q75), color = "#FAC898", size = 2) + 
  geom_point(fill = "#FAC898", size = 5, shape = 21, color = "black") +
  labs(y = "", x = expression("Benthivory ("*gC~m^{-2}~d^{-1}*")")) + theme_minimal(base_size = 14) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_S4F <- ridge_summary %>% dplyr::filter(Function == "Multifunctionality") |> 
  ggplot(aes(y = fct_rev(factor(period)), x = mean)) +
  geom_linerange(aes(xmin = q05, xmax = q95), color = "#f4cccc", size = 1) +   scale_x_log10() +
  geom_linerange(aes(xmin = q25, xmax = q75), color = "#FF968A", size = 2) + 
  geom_point(fill = "#FF968A", size = 5, shape = 21, color = "black") +
  labs(y = "", x = "Multifunctional Index") + theme_minimal(base_size = 14) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

# Altogether
Figure_S4 = Figure_S4A + Figure_S4B + Figure_S4C + Figure_S4D + Figure_S4E + Figure_S4F + 
  plot_layout(guides = "collect", ncol = 2) 

#### Export the data  ----
## Figures
ggsave(Figure_S3, filename = "Figure_S3.png", path = "Outputs/", device = "png", width = 12, height = 9, dpi = 300)  
ggsave(Figure_S4, filename = "Figure_S4.png", path = "Outputs/", device = "png", width = 10, height = 12, dpi = 300)  