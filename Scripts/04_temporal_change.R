#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
Hexagonal_grid     <- st_read("Data/Grid_0-1000m_Med.shp")
Guilds             <- readxl::read_xlsx("Data/Guilds_MED.xlsx") |> mutate(SPECIES = str_replace_all(SPECIES, "\u00A0", " "))

## Charge from previous scripts
load("Outputs/dat_proc/Med_all.RData")
load("Outputs/dat_proc/Medit_Temp.RData")
load("Outputs/dat_proc/medits_sf_percentile.Rdata")
load("Outputs/dat_proc/Medit_FunCatch_without_NA.RData")

## Color palette and map
land       <- ne_countries(scale = "medium", returnclass = "sf")

## Functions
source("Scripts/00_functions_script.R")

## Check visually hexagon grid with data
Hexagonal_grid_with_points <- st_join(st_transform(Hexagonal_grid, st_crs(medits_sf_percentile)), 
                                      medits_sf_percentile, join = st_intersects) |> drop_na(YEAR) |> 
  ggplot() + geom_sf(fill = "red", alpha = 0.4, color = "black") + theme_minimal() +
  geom_sf(data = ne_countries(scale = "medium", returnclass = "sf"), fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

## Attribute an ID from Hexagonal_grid to dataset
medits_with_hexID <- st_join(medits_sf_percentile, st_transform(Hexagonal_grid, st_crs(medits_sf_percentile)), 
  join = st_intersects) |> drop_na(id) |> rename(HEX_ID = id) |> relocate(HEX_ID, .after = YEAR)

## Keep only HEX cells with at least 10 years of data
medits_with_hexID_over_10years <- medits_with_hexID |> 
  filter(HEX_ID %in% (medits_with_hexID |> group_by(HEX_ID) |> summarise(n_years = n_distinct(YEAR)) |> 
                        filter(n_years >= 10) |> st_drop_geometry() |> ungroup() |> select(HEX_ID) |> 
                        pull(HEX_ID))) |> group_by(HEX_ID) |> 
  dplyr::select(c(1:3, 5:9, 15, 22:27, 29:32))

## Analyze slopes over year
Slopes <- st_as_sf(medits_with_hexID_over_10years |>
  pivot_longer(cols = c(community_Gc, community_Fn, community_Fp, Ic_plank, Ic_benthivorous, Multifunctionality),
    names_to = "metric", values_to = "value") |> group_by(HEX_ID, metric) |> do(tidy(lm(value ~ YEAR, data = .))) |>
  filter(term == "YEAR") |> select(HEX_ID, metric, slope = estimate, p.value = p.value) |>
  pivot_wider(names_from = metric, values_from = c(slope, p.value), names_glue = "{metric}_{.value}") |>
  left_join(medits_with_hexID_over_10years |> st_drop_geometry() |> group_by(HEX_ID) |>
      summarise(centroid_lon = mean((left + right) / 2, na.rm = TRUE),
                centroid_lat = mean((top + bottom) / 2, na.rm = TRUE)) |>
      st_as_sf(coords = c("centroid_lon", "centroid_lat"), crs = st_crs(medits_with_hexID_over_10years))))  |> 
  mutate(Fn_trend_class  = case_when(community_Fn_slope       < 0 & community_Fn_p.value       <= 0.05 ~ "1. significantly negative",
                                     community_Fn_slope       < 0 & community_Fn_p.value       >  0.05 ~ "2. negative trend",
                                     community_Fn_slope       > 0 & community_Fn_p.value       >  0.05 ~ "3. positive trend",
                                     community_Fn_slope       > 0 & community_Fn_p.value       <= 0.05 ~ "4. significantly positive"),
         Fp_trend_class  = case_when(community_Fp_slope       < 0 & community_Fp_p.value       <= 0.05 ~ "1. significantly negative",
                                     community_Fp_slope       < 0 & community_Fp_p.value       >  0.05 ~ "2. negative trend",
                                     community_Fp_slope       > 0 & community_Fp_p.value       >  0.05 ~ "3. positive trend",
                                     community_Fp_slope       > 0 & community_Fp_p.value       <= 0.05 ~ "4. significantly positive"),
         Gc_trend_class  = case_when(community_Gc_slope       < 0 & community_Gc_p.value       <= 0.05 ~ "1. significantly negative",
                                     community_Gc_slope       < 0 & community_Gc_p.value       >  0.05 ~ "2. negative trend",
                                     community_Gc_slope       > 0 & community_Gc_p.value       >  0.05 ~ "3. positive trend",
                                     community_Gc_slope       > 0 & community_Gc_p.value       <= 0.05 ~ "4. significantly positive"),
         IcB_trend_class = case_when(Ic_benthivorous_slope    < 0 & Ic_benthivorous_p.value    <= 0.05 ~ "1. significantly negative",
                                     Ic_benthivorous_slope    < 0 & Ic_benthivorous_p.value    >  0.05 ~ "2. negative trend",
                                     Ic_benthivorous_slope    > 0 & Ic_benthivorous_p.value    >  0.05 ~ "3. positive trend",
                                     Ic_benthivorous_slope    > 0 & Ic_benthivorous_p.value    <= 0.05 ~ "4. significantly positive"),
         IcP_trend_class = case_when(Ic_plank_slope           < 0 & Ic_plank_p.value           <= 0.05 ~ "1. significantly negative",
                                     Ic_plank_slope           < 0 & Ic_plank_p.value           >  0.05 ~ "2. negative trend",
                                     Ic_plank_slope           > 0 & Ic_plank_p.value           >  0.05 ~ "3. positive trend",
                                     Ic_plank_slope           > 0 & Ic_plank_p.value           <= 0.05 ~ "4. significantly positive"),
         Mf_trend_class  = case_when(Multifunctionality_slope < 0 & Multifunctionality_p.value <= 0.05 ~ "1. significantly negative",
                                     Multifunctionality_slope < 0 & Multifunctionality_p.value >  0.05 ~ "2. negative trend",
                                     Multifunctionality_slope > 0 & Multifunctionality_p.value >  0.05 ~ "3. positive trend",
                                     Multifunctionality_slope > 0 & Multifunctionality_p.value <= 0.05 ~ "4. significantly positive")) 

## Plots
Figure_4A1 <- ggplot() +
  geom_sf(data = subset(Slopes, Gc_trend_class %in% c("2. negative trend", "3. positive trend")),
          aes(fill = Gc_trend_class, color = Gc_trend_class, size = Gc_trend_class, shape = Gc_trend_class)) +
  geom_sf(data = subset(Slopes, Gc_trend_class %in% c("1. significantly negative", "4. significantly positive")),
          aes(fill = Gc_trend_class, color = Gc_trend_class, size = Gc_trend_class, shape = Gc_trend_class)) +
  scale_shape_manual(values = c("1. significantly negative" = 21, "4. significantly positive" = 21, 
                                "2. negative trend" = 20, "3. positive trend" = 20)) +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "4. significantly positive" = "#9999FF", 
                               "2. negative trend" = "#000000", "3. positive trend" = "#000000")) +
  scale_color_manual(values = c("1. significantly negative" = "#000000", "4. significantly positive" = "#000000", 
                                "2. negative trend" = "#FF9999", "3. positive trend" = "#CCCCFF")) +
  scale_size_manual(values = c("1. significantly negative" = 4, "4. significantly positive" = 4, 
                               "2. negative trend" = 2, "3. positive trend" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Production")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_4B1 <- ggplot() +
  geom_sf(data = subset(Slopes, Fn_trend_class %in% c("2. negative trend", "3. positive trend")),
          aes(fill = Fn_trend_class, color = Fn_trend_class, size = Fn_trend_class, shape = Fn_trend_class)) +
  geom_sf(data = subset(Slopes, Fn_trend_class %in% c("1. significantly negative", "4. significantly positive")),
          aes(fill = Fn_trend_class, color = Fn_trend_class, size = Fn_trend_class, shape = Fn_trend_class)) +
  scale_shape_manual(values = c("1. significantly negative" = 21, "4. significantly positive" = 21, 
                                "2. negative trend" = 20, "3. positive trend" = 20)) +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "4. significantly positive" = "#9999FF", 
                               "2. negative trend" = "#000000", "3. positive trend" = "#000000")) +
  scale_color_manual(values = c("1. significantly negative" = "#000000", "4. significantly positive" = "#000000", 
                                "2. negative trend" = "#FF9999", "3. positive trend" = "#CCCCFF")) +
  scale_size_manual(values = c("1. significantly negative" = 4, "4. significantly positive" = 4, 
                               "2. negative trend" = 2, "3. positive trend" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Nitrogen excretion")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_4C1 <- ggplot() +
  geom_sf(data = subset(Slopes, Fp_trend_class %in% c("2. negative trend", "3. positive trend")),
          aes(fill = Fp_trend_class, color = Fp_trend_class, size = Fp_trend_class, shape = Fp_trend_class)) +
  geom_sf(data = subset(Slopes, Fp_trend_class %in% c("1. significantly negative", "4. significantly positive")),
          aes(fill = Fp_trend_class, color = Fp_trend_class, size = Fp_trend_class, shape = Fp_trend_class)) +
  scale_shape_manual(values = c("1. significantly negative" = 21, "4. significantly positive" = 21, 
                                "2. negative trend" = 20, "3. positive trend" = 20)) +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "4. significantly positive" = "#9999FF", 
                               "2. negative trend" = "#000000", "3. positive trend" = "#000000")) +
  scale_color_manual(values = c("1. significantly negative" = "#000000", "4. significantly positive" = "#000000", 
                                "2. negative trend" = "#FF9999", "3. positive trend" = "#CCCCFF")) +
  scale_size_manual(values = c("1. significantly negative" = 4, "4. significantly positive" = 4, 
                               "2. negative trend" = 2, "3. positive trend" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Phosphorus excretion")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_4D1 <- ggplot() +
  geom_sf(data = subset(Slopes, IcP_trend_class %in% c("2. negative trend", "3. positive trend")),
          aes(fill = IcP_trend_class, color = IcP_trend_class, size = IcP_trend_class, shape = IcP_trend_class)) +
  geom_sf(data = subset(Slopes, IcP_trend_class %in% c("1. significantly negative", "4. significantly positive")),
          aes(fill = IcP_trend_class, color = IcP_trend_class, size = IcP_trend_class, shape = IcP_trend_class)) +
  scale_shape_manual(values = c("1. significantly negative" = 21, "4. significantly positive" = 21, 
                                "2. negative trend" = 20, "3. positive trend" = 20)) +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "4. significantly positive" = "#9999FF", 
                               "2. negative trend" = "#000000", "3. positive trend" = "#000000")) +
  scale_color_manual(values = c("1. significantly negative" = "#000000", "4. significantly positive" = "#000000", 
                                "2. negative trend" = "#FF9999", "3. positive trend" = "#CCCCFF")) +
  scale_size_manual(values = c("1. significantly negative" = 4, "4. significantly positive" = 4, 
                               "2. negative trend" = 2, "3. positive trend" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Planktivory")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_4E1 <- ggplot() +
  geom_sf(data = subset(Slopes, IcB_trend_class %in% c("2. negative trend", "3. positive trend")),
          aes(fill = IcB_trend_class, color = IcB_trend_class, size = IcB_trend_class, shape = IcB_trend_class)) +
  geom_sf(data = subset(Slopes, IcB_trend_class %in% c("1. significantly negative", "4. significantly positive")),
          aes(fill = IcB_trend_class, color = IcB_trend_class, size = IcB_trend_class, shape = IcB_trend_class)) +
  scale_shape_manual(values = c("1. significantly negative" = 21, "4. significantly positive" = 21, 
                                "2. negative trend" = 20, "3. positive trend" = 20)) +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "4. significantly positive" = "#9999FF", 
                               "2. negative trend" = "#000000", "3. positive trend" = "#000000")) +
  scale_color_manual(values = c("1. significantly negative" = "#000000", "4. significantly positive" = "#000000", 
                                "2. negative trend" = "#FF9999", "3. positive trend" = "#CCCCFF")) +
  scale_size_manual(values = c("1. significantly negative" = 4, "4. significantly positive" = 4, 
                               "2. negative trend" = 2, "3. positive trend" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Benthivory")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_4F1 <- ggplot() +
  geom_sf(data = subset(Slopes, Mf_trend_class %in% c("2. negative trend", "3. positive trend")),
          aes(fill = Mf_trend_class, color = Mf_trend_class, size = Mf_trend_class, shape = Mf_trend_class)) +
  geom_sf(data = subset(Slopes, Mf_trend_class %in% c("1. significantly negative", "4. significantly positive")),
          aes(fill = Mf_trend_class, color = Mf_trend_class, size = Mf_trend_class, shape = Mf_trend_class)) +
  scale_shape_manual(values = c("1. significantly negative" = 21, "4. significantly positive" = 21, 
                                "2. negative trend" = 20, "3. positive trend" = 20)) +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "4. significantly positive" = "#9999FF", 
                               "2. negative trend" = "#000000", "3. positive trend" = "#000000")) +
  scale_color_manual(values = c("1. significantly negative" = "#000000", "4. significantly positive" = "#000000", 
                                "2. negative trend" = "#FF9999", "3. positive trend" = "#CCCCFF")) +
  scale_size_manual(values = c("1. significantly negative" = 4, "4. significantly positive" = 4, 
                               "2. negative trend" = 2, "3. positive trend" = 2)) +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  theme_minimal() + coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) + 
  ggtitle(expression("Multifunctional Index")) +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

Figure_4A2 <- Slopes %>% st_drop_geometry() %>% group_by(Class = Gc_trend_class) %>% summarise(n = n()) |> 
  mutate(Fraction = n / sum(n)) |> ggplot(aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "2. negative trend" = "#FF9999",
                               "3. positive trend" = "#CCCCFF", "4. significantly positive" = "#9999FF")) +
  theme(legend.position = "none")

Figure_4B2 <- Slopes %>% st_drop_geometry() %>% group_by(Class = Fn_trend_class) %>% summarise(n = n()) |> 
  mutate(Fraction = n / sum(n)) |> ggplot(aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "2. negative trend" = "#FF9999",
                               "3. positive trend" = "#CCCCFF", "4. significantly positive" = "#9999FF")) +
  theme(legend.position = "none")

Figure_4C2 <- Slopes %>% st_drop_geometry() %>% group_by(Class = Fp_trend_class) %>% summarise(n = n()) |> 
  mutate(Fraction = n / sum(n)) |> ggplot(aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "2. negative trend" = "#FF9999",
                               "3. positive trend" = "#CCCCFF", "4. significantly positive" = "#9999FF")) +
  theme(legend.position = "none")

Figure_4D2 <- Slopes %>% st_drop_geometry() %>% group_by(Class = IcP_trend_class) %>% summarise(n = n()) |> 
  mutate(Fraction = n / sum(n)) |> ggplot(aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "2. negative trend" = "#FF9999",
                               "3. positive trend" = "#CCCCFF", "4. significantly positive" = "#9999FF")) +
  theme(legend.position = "none")

Figure_4E2 <- Slopes %>% st_drop_geometry() %>% group_by(Class = IcB_trend_class) %>% summarise(n = n()) |> 
  mutate(Fraction = n / sum(n)) |> ggplot(aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "2. negative trend" = "#FF9999",
                               "3. positive trend" = "#CCCCFF", "4. significantly positive" = "#9999FF")) +
  theme(legend.position = "none")

Figure_4F2 <- Slopes %>% st_drop_geometry() %>% group_by(Class = Mf_trend_class) %>% summarise(n = n()) |> 
  mutate(Fraction = n / sum(n)) |> ggplot(aes(x = 2, y = Fraction, fill = Class)) +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + theme_void() +
  scale_fill_manual(values = c("1. significantly negative" = "#FF6666", "2. negative trend" = "#FF9999",
                               "3. positive trend" = "#CCCCFF", "4. significantly positive" = "#9999FF")) +
  theme(legend.position = "none")

  # Add pie plot insets
Figure_4A <- Figure_4A1 + annotation_custom(grob = ggplotGrob(Figure_4A2), xmin = 25, xmax = 40, ymin = 40, ymax = 46)
Figure_4B <- Figure_4B1 + annotation_custom(grob = ggplotGrob(Figure_4B2), xmin = 25, xmax = 40, ymin = 40, ymax = 46)
Figure_4C <- Figure_4C1 + annotation_custom(grob = ggplotGrob(Figure_4C2), xmin = 25, xmax = 40, ymin = 40, ymax = 46)
Figure_4D <- Figure_4D1 + annotation_custom(grob = ggplotGrob(Figure_4D2), xmin = 25, xmax = 40, ymin = 40, ymax = 46)
Figure_4E <- Figure_4E1 + annotation_custom(grob = ggplotGrob(Figure_4E2), xmin = 25, xmax = 40, ymin = 40, ymax = 46)
Figure_4F <- Figure_4F1 + annotation_custom(grob = ggplotGrob(Figure_4F2), xmin = 25, xmax = 40, ymin = 40, ymax = 46)

# Final Figure
Figure_4 = Figure_4A + Figure_4B + Figure_4C + Figure_4D + Figure_4E + Figure_4F + 
  plot_layout(guides = "collect", ncol = 2) 

#### Export the data  ----
## Figures
ggsave(Figure_4, filename = "Figure_4.png", path = "Outputs/", device = "png", width = 12, height = 9, dpi = 300)  