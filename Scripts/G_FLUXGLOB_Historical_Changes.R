#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("rphylopic") ; library("FNN") ; library("pbapply")

## Charge from previous scripts
Feeding <- readxl::read_excel("Data/Feeding_mode_filled.xlsx")
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB_with_fluxes.RData")

## Functions
source("Scripts/00_functions_script.R")

## Merge Final dataset with feeding guilds
## One species has empty space
FLUXGLOB_final$Species_FLUXGLOB[FLUXGLOB_final$Species_FLUXGLOB == "Symphodus doderleini "] = "Symphodus doderleini"
FLUXGLOB_final = FLUXGLOB_final |> left_join(Feeding)

# Community computation at the individual level
FLUXGLOB_community = FLUXGLOB_final |> 
  mutate(Biomass = Weight * Num / Haul_swept_area / (Haul_duration / 60 / 24),
         community_Fn = Fn_mean * Num / Haul_swept_area / (Haul_duration / 60 / 24),
         community_Fp = Fp_mean * Num / Haul_swept_area / (Haul_duration / 60 / 24),
         community_Gc = Gc_mean * Num / Haul_swept_area / (Haul_duration / 60 / 24)) |> 
  group_by(Haul_ID, cell, Depth, SBTemp, ID_Fig1, Year, Month, geometry) |> 
  summarise(Biomass = sum(Biomass), community_Fn = sum(community_Fn), community_Fp = sum(community_Fp),
            community_Gc = sum(community_Gc), .groups = "drop") 

# Community computation at the trophic guild level
FLUXGLOB_community_comm = FLUXGLOB_final |> 
  group_by(Haul_ID, cell, Depth, SBTemp, ID_Fig1, Year, Month, geometry, Feeding_mode) |> 
  summarise(community_Ic = sum(Ic_mean), .groups = "drop") |> 
  dplyr::filter(Feeding_mode %in% c("Benthivorous", "Generalist", "Planktivorous", "Piscivorous")) |> 
  mutate(Ic_plank = case_when(Feeding_mode == "Planktivorous" ~ community_Ic,
                              Feeding_mode == "Generalist" ~ 0.25 * community_Ic, TRUE ~ 0),
         Ic_benth = case_when(Feeding_mode == "Benthivorous"  ~ community_Ic,
                              Feeding_mode == "Generalist" ~ 0.25 * community_Ic, TRUE ~ 0),
         Ic_pisci = case_when(Feeding_mode == "Piscivorous"  ~ community_Ic,
                              Feeding_mode == "Generalist" ~ 0.25 * community_Ic, TRUE ~ 0)) |> 
  group_by(Haul_ID, cell, Depth, SBTemp, ID_Fig1, Year, Month, geometry) |> 
  summarise(Ic_plank = sum(Ic_plank), Ic_benth = sum(Ic_benth), Ic_pisci = sum(Ic_pisci), .groups = "drop") 

# Merge all results
FLUXGLOB_community <- FLUXGLOB_community |> left_join(FLUXGLOB_community_comm |> sf::st_drop_geometry())
rm(FLUXGLOB_community_comm)

# Let's plot the last years
last_years <- FLUXGLOB_community |> st_drop_geometry() |> group_by(ID_Fig1) |>
  summarise(last_year = max(Year), .groups = "drop")
FLUX_years <- FLUXGLOB_community |> left_join(last_years, by = "ID_Fig1") |>
  mutate(Period = if_else(Year >= last_year - 4, "Recent", "Historical"))

FLUX_year <- FLUX_years |> st_drop_geometry() |> group_by(ID_Fig1, Year, cell, Period) |>
  summarise(Fn_geo    = exp(mean(log(community_Fn + 1), na.rm = TRUE)) - 1,
            Fp_geo    = exp(mean(log(community_Fp + 1), na.rm = TRUE)) - 1,
            Gc_geo    = exp(mean(log(community_Gc + 1), na.rm = TRUE)) - 1,
            IcPl_geo  = exp(mean(log(Ic_plank + 1), na.rm = TRUE)) - 1,
            IcPi_geo  = exp(mean(log(Ic_pisci + 1), na.rm = TRUE)) - 1,
            IcB_geo   = exp(mean(log(Ic_benth + 1), na.rm = TRUE)) - 1,
            n_hauls   = n(), .groups = "drop") |> 
  mutate(across(all_of(c("Fn_geo", "Fp_geo", "Gc_geo", "IcPl_geo", "IcPi_geo", "IcB_geo")),
      ~ pmin(.x / quantile(.x, 0.99, na.rm = TRUE), 1), .names = "{.col}_norm")) |> rowwise() |> 
  mutate(Multifunctionality = exp(mean(log(c_across(ends_with("_norm")) + 1e-6), na.rm = TRUE))) |> ungroup()

Flux_long        <- FLUX_year |> dplyr::select(ID_Fig1, Year, cell, Period, Fn_geo, Fp_geo, Gc_geo, IcPl_geo, IcPi_geo, IcB_geo, Multifunctionality) |>
  pivot_longer(cols = c(Fn_geo, Fp_geo, Gc_geo, IcPl_geo, IcPi_geo, IcB_geo, Multifunctionality), names_to = "Flux", values_to = "Value")
Baseline         <- Flux_long |> filter(Period == "Historical") |> group_by(ID_Fig1, Flux) |>
  summarise(Mean_hist = mean(Value, na.rm = TRUE), SD_hist = sd(Value, na.rm = TRUE), .groups = "drop")
Flux_z           <- Flux_long |> left_join(Baseline, by = c("ID_Fig1","Flux")) |> mutate(Z = (Value - Mean_hist) / SD_hist) |>
  mutate(ID_plot = factor(as.character(ID_Fig1), levels = rev(as.character(1:38))))
Recent_summary   <- Flux_z |> filter(Period == "Recent") |> group_by(ID_Fig1, Flux) |>
  summarise(median_recent = median(Z, na.rm = TRUE), q25_recent = quantile(Z, 0.25, na.rm = TRUE),
            q75_recent    = quantile(Z, 0.75, na.rm = TRUE), .groups = "drop") |>
  mutate(ID_plot = factor(as.character(ID_Fig1), levels = rev(as.character(1:38)))) |> 
  mutate(across(everything(), ~replace_na(., 0))) |> 
  mutate(significance = factor(case_when(q25_recent < 0 & q75_recent < 0 ~ "Negative", q25_recent > 0 & 
                                           q75_recent > 0 ~ "Positive", TRUE ~ "Not significant"),
                               levels = c("Negative", "Not significant", "Positive")))

(Figure_3A = ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
  geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
  geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
  geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
  geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
  geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
  geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Gc_geo"), 
               aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
  geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "Gc_geo"), Period == "Historical"), aes(x = Z, y = ID_plot),
             colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
  geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Gc_geo"), 
               aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
               color = "#0047AB", linewidth = 1) +
  geom_point(data = Recent_summary |> dplyr::filter(Flux == "Gc_geo"), 
             aes(x = median_recent, y = ID_plot, fill = significance), colour = "#0047AB", size = 3.5, 
             color = "black", shape = 21, show.legend = F) +
  scale_fill_manual(values = c("#0047AB", "white", "#6495ED")) +
  scale_y_discrete(breaks = rev(as.character(1:38)),
                   labels = function(x) ifelse(as.numeric(as.character(x)) %% 3 == 1, x, "")) +
  labs(x = "", y = "Ecoregion", title = expression(bold("A.") ~ "Carbon Production")) +
  theme_classic(base_size = 13) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.8, "cm"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

(Figure_3B = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
    geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Fn_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
    geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "Fn_geo"), Period == "Historical"), aes(x = Z, y = ID_plot),
               colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Fn_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
                 color = "#4F7942", linewidth = 1) +
    geom_point(data = Recent_summary |> dplyr::filter(Flux == "Fn_geo"), 
               aes(x = median_recent, y = ID_plot, fill = significance), colour = "#4F7942", size = 3.5, 
               color = "black", shape = 21, show.legend = F) +
    scale_fill_manual(values = c("#4F7942", "white", "#93C572")) +
    labs(x = "", y = "Ecoregion", title = expression(bold("B.") ~ "Nitrogen excretion")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "cm"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

(Figure_3C = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
    geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Fp_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
    geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "Fp_geo"), Period == "Historical"), aes(x = Z, y = ID_plot),
               colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Fp_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
                 color = "#F4C430", linewidth = 1) +
    geom_point(data = Recent_summary |> dplyr::filter(Flux == "Fp_geo"), 
               aes(x = median_recent, y = ID_plot, fill = significance), colour = "#F4C430", size = 3.5, 
               color = "black", shape = 21, show.legend = F) +
    scale_fill_manual(values = c("#F4C430", "white", "#FBEC5D")) +
    labs(x = "Standardized flux anomaly (σ)", y = "Ecoregion", title = expression(bold("C.") ~ "Phosphorus excretion")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "cm"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

(Figure_3D = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
    geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "IcB_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
    geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "IcB_geo"), Period == "Historical"), aes(x = Z, y = ID_plot),
               colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "IcB_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
                 color = "#A300A3", linewidth = 1) +
    geom_point(data = Recent_summary |> dplyr::filter(Flux == "IcB_geo"), 
               aes(x = median_recent, y = ID_plot, fill = significance), colour = "#A300A3", size = 3.5, 
               color = "black", shape = 21, show.legend = F) +
    scale_fill_manual(values = c("#A300A3", "white", "#CF9FFF")) +
    scale_y_discrete(breaks = rev(as.character(1:38)),
                     labels = function(x) ifelse(as.numeric(as.character(x)) %% 3 == 1, x, "")) +
    labs(x = "", y = "Ecoregion", title = expression(bold("D.") ~ "Benthivory")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "cm"),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 13),
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 15),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

(Figure_3E = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
    geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "IcPi_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
    geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "IcPi_geo"), Period == "Historical"), aes(x = Z, y = ID_plot),
               colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "IcPi_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
                 color = "#708090", linewidth = 1) +
    geom_point(data = Recent_summary |> dplyr::filter(Flux == "IcPi_geo"), 
               aes(x = median_recent, y = ID_plot, fill = significance), colour = "#708090", size = 3.5, 
               color = "black", shape = 21, show.legend = F) +
    scale_fill_manual(values = c("#708090", "white", "#7393B3")) +
    labs(x = "Standardized flux anomaly (σ)", y = "Ecoregion", title = expression(bold("E.") ~ "Piscivory")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "cm"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 15),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

(Figure_3F = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
    geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "IcPl_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
    geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "IcPl_geo"), Period == "Historical"), aes(x = Z, y = ID_plot),
               colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "IcPl_geo"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
                 color = "#D27D2D", linewidth = 1) +
    geom_point(data = Recent_summary |> dplyr::filter(Flux == "IcPl_geo"), 
               aes(x = median_recent, y = ID_plot, fill = significance), colour = "#D27D2D", size = 3.5, 
               color = "black", shape = 21, show.legend = F) +
    scale_fill_manual(values = c("#D27D2D", "white", "#FFAC1C")) +
    labs(x = "", y = "Ecoregion", title = expression(bold("F.") ~ "Planktivory")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "cm"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 15),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

(Figure_3G = ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + coord_cartesian(xlim = c(-3, 3)) +
    geom_hline(yintercept = 7.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 9.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 12.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 16.5, linetype = "dotted", colour = "black") +
    geom_hline(yintercept = 29.5, linetype = "dotted", colour = "black") +    
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Multifunctionality"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot), linewidth = 0) +
    geom_point(data = filter(Flux_z |> dplyr::filter(Flux == "Multifunctionality"), Period == "Historical"), aes(x = Z, y = ID_plot),
               colour = "grey40", alpha = 0.05, size = 0.9, position = position_jitter(width = 0, height = 0.12, seed = 123)) +
    geom_segment(data = Recent_summary |> dplyr::filter(Flux == "Multifunctionality"), 
                 aes(x = q25_recent, xend = q75_recent, y = ID_plot, yend = ID_plot),
                 color = "#DC143C", linewidth = 1) +
    geom_point(data = Recent_summary |> dplyr::filter(Flux == "Multifunctionality"), 
               aes(x = median_recent, y = ID_plot, fill = significance), colour = "#DC143C", size = 3.5, 
               color = "black", shape = 21, show.legend = F) +
    scale_fill_manual(values = c("#DC143C", "white", "#F88379")) +
    scale_y_discrete(breaks = rev(as.character(1:38)),
                     labels = function(x) ifelse(as.numeric(as.character(x)) %% 2 == 1, x, "")) +
    labs(x = "Standardized anomaly (σ)", y = "Ecoregion", title = expression(bold("G.") ~ "Multi-Fluxes Index")) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "cm"),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_text(size = 13),
          axis.text.x = element_text(size = 13),
          axis.title.x = element_text(size = 15),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)))

design <- "ABCGG
           DEFGG"

(Figure_3 <- Figure_3A + Figure_3B + Figure_3C + Figure_3D + Figure_3E + Figure_3F + Figure_3G + plot_layout(design = design))

ggsave(Figure_3, filename = "FLUXGLOB/Raw/Figure_3_Eli.png", path = "Outputs/", device = "png", width = 14, height = 8, dpi = 300)  
