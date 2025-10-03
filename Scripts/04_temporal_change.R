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
my_palette <- c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffff",
                "#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")

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
                        pull(HEX_ID))) |> group_by(HEX_ID) 

## Analyze slopes over year
Slopes <- st_as_sf(medits_with_hexID_over_10years |>
  pivot_longer(cols = c(community_Gc, community_Fn, community_Fp, Ic_plank, Ic_benthivorous),
    names_to = "metric", values_to = "value") |> group_by(HEX_ID, metric) |> do(tidy(lm(value ~ YEAR, data = .))) |>
  filter(term == "YEAR") |> select(HEX_ID, metric, slope = estimate) |>
  pivot_wider(names_from = metric, values_from = slope, names_glue = "{metric}_slope") |>
  left_join(medits_with_hexID_over_10years |> st_drop_geometry() |> group_by(HEX_ID) |>
      summarise(centroid_lon = mean((left + right) / 2, na.rm = TRUE),
                centroid_lat = mean((top + bottom) / 2, na.rm = TRUE)) |>
      st_as_sf(coords = c("centroid_lon", "centroid_lat"), crs = st_crs(medits_with_hexID_over_10years)))) 

## Define class
# Carbon
qs <- quantile(Slopes$community_Gc_slope, probs = c(0.01, 0.99), na.rm = TRUE) ; q1  <- qs[1]; q99 <- qs[2]
Slopes <- Slopes |> mutate(Gc_trend_class = case_when(community_Gc_slope < q1  ~ "<1%",
  community_Gc_slope >= q1 & community_Gc_slope <= q99 & community_Gc_slope < 0 ~ "Negative",
  community_Gc_slope >= q1 & community_Gc_slope <= q99 & community_Gc_slope > 0 ~ "Positive",
  community_Gc_slope > q99 ~ ">99%")) 
# Nitrogen
qs <- quantile(Slopes$community_Fn_slope, probs = c(0.01, 0.99), na.rm = TRUE) ; q1  <- qs[1]; q99 <- qs[2]
Slopes <- Slopes |> mutate(Fn_trend_class = case_when(community_Fn_slope < q1  ~ "<1%",
  community_Fn_slope >= q1 & community_Fn_slope <= q99 & community_Fn_slope < 0 ~ "Negative",
  community_Fn_slope >= q1 & community_Fn_slope <= q99 & community_Fn_slope > 0 ~ "Positive",
  community_Fn_slope > q99 ~ ">99%"))
# Phosphorus
qs <- quantile(Slopes$community_Fp_slope, probs = c(0.01, 0.99), na.rm = TRUE) ; q1  <- qs[1]; q99 <- qs[2]
Slopes <- Slopes |> mutate(Fp_trend_class = case_when(community_Fp_slope < q1  ~ "<1%",
  community_Fp_slope >= q1 & community_Fp_slope <= q99 & community_Fp_slope < 0 ~ "Negative",
  community_Fp_slope >= q1 & community_Fp_slope <= q99 & community_Fp_slope > 0 ~ "Positive",
  community_Fp_slope > q99 ~ ">99%"))
# Planktivory
qs <- quantile(Slopes$Ic_plank_slope, probs = c(0.01, 0.99), na.rm = TRUE) ; q1  <- qs[1]; q99 <- qs[2]
Slopes <- Slopes |> mutate(IcP_trend_class = case_when(Ic_plank_slope < q1  ~ "<1%",
  Ic_plank_slope >= q1 & Ic_plank_slope <= q99 & Ic_plank_slope < 0 ~ "Negative",
  Ic_plank_slope >= q1 & Ic_plank_slope <= q99 & Ic_plank_slope > 0 ~ "Positive",
  Ic_plank_slope > q99 ~ ">99%"))
# Benthivory
qs <- quantile(Slopes$Ic_benthivorous_slope, probs = c(0.01, 0.99), na.rm = TRUE) ; q1  <- qs[1]; q99 <- qs[2]
Slopes <- Slopes |> mutate(IcB_trend_class = case_when(Ic_benthivorous_slope < q1  ~ "<1%",
  Ic_benthivorous_slope >= q1 & Ic_benthivorous_slope <= q99 & Ic_benthivorous_slope < 0 ~ "Negative",
  Ic_benthivorous_slope >= q1 & Ic_benthivorous_slope <= q99 & Ic_benthivorous_slope > 0 ~ "Positive",
  Ic_benthivorous_slope > q99 ~ ">99%"))

# Reorder the slopes
Slopes <- Slopes |> mutate(
    Fn_trend_class = factor(Fn_trend_class, levels = c("<1%", "Negative", "Positive", ">99%")),
    Gc_trend_class = factor(Gc_trend_class, levels = c("<1%", "Negative", "Positive", ">99%")),
    Fp_trend_class = factor(Fp_trend_class, levels = c("<1%", "Negative", "Positive", ">99%")),
    IcP_trend_class = factor(IcP_trend_class, levels = c("<1%", "Negative", "Positive", ">99%")),
    IcB_trend_class = factor(IcB_trend_class, levels = c("<1%", "Negative", "Positive", ">99%")))

Figure_3A <- ggplot(Slopes) + geom_sf(aes(fill = Gc_trend_class), color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("<1%" = "#d73027", "Negative" = "#E88B95", 
                               "Positive" = "#abd9e9", ">99%" = "#4169e1"), name = "") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme_minimal() + ggtitle("Fish Production") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_3B <- ggplot(Slopes) + geom_sf(aes(fill = Fn_trend_class), color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("<1%" = "#d73027", "Negative" = "#E88B95", 
                               "Positive" = "#abd9e9", ">99%" = "#4169e1"), name = "") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme_minimal() + ggtitle("Fish Nitrogen Excretion") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_3C <- ggplot(Slopes) + geom_sf(aes(fill = Fp_trend_class), color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("<1%" = "#d73027", "Negative" = "#E88B95", 
                               "Positive" = "#abd9e9", ">99%" = "#4169e1"), name = "") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme_minimal() + ggtitle("Fish Phosphorus Excretion") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_3D <- ggplot(Slopes) + geom_sf(aes(fill = IcP_trend_class), color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("<1%" = "#d73027", "Negative" = "#E88B95", 
                               "Positive" = "#abd9e9", ">99%" = "#4169e1"), name = "") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme_minimal() + ggtitle("Fish Planktivory") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

Figure_3E <- ggplot(Slopes) + geom_sf(aes(fill = IcB_trend_class), color = "black", shape = 21, size = 4) +
  scale_fill_manual(values = c("<1%" = "#d73027", "Negative" = "#E88B95", 
                               "Positive" = "#abd9e9", ">99%" = "#4169e1"), name = "") +
  geom_sf(data = land, fill = "lightgray", color = "black") +
  coord_sf(xlim = c(-5, 35), ylim = c(34, 46)) +
  theme_minimal() + ggtitle("Fish Benthivory") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12),
        legend.position = "bottom")

# Classify slopes into Positive / Negative
doughnut_data <- Slopes %>% st_drop_geometry() %>%
  select(all_of(c("community_Gc_slope", "community_Fn_slope", "community_Fp_slope", 
                          "Ic_plank_slope", "Ic_benthivorous_slope"))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Slope") %>%
  mutate(Sign = case_when(Slope > 0 ~ "Positive", Slope < 0 ~ "Negative", TRUE ~ "Zero")) %>%
  filter(Sign != "Zero") %>% group_by(Variable, Sign) %>% summarise(Count = n(), .groups = "drop") %>%
  group_by(Variable) %>% mutate(Fraction = Count / sum(Count))

# Doughnut plot
Gc_doughnut <- doughnut_data |> filter(Variable == "community_Gc_slope") |> 
  ggplot(aes(x = 2, y = Fraction, fill = Sign)) + ggtitle("") +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + 
  scale_fill_manual(values = c("Positive" = "#abd9e9", "Negative" = "#E88B95")) +
  theme_void() + theme(legend.position = "none", strip.text = element_text(size = 12))
Fn_doughnut <- doughnut_data |> filter(Variable == "community_Fn_slope") |> 
  ggplot(aes(x = 2, y = Fraction, fill = Sign)) + ggtitle("") +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + 
  scale_fill_manual(values = c("Positive" = "#abd9e9", "Negative" = "#E88B95")) +
  theme_void() + theme(legend.position = "none", strip.text = element_text(size = 12))
Fp_doughnut <- doughnut_data |> filter(Variable == "community_Fp_slope") |> 
  ggplot(aes(x = 2, y = Fraction, fill = Sign)) + ggtitle("") +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + 
  scale_fill_manual(values = c("Positive" = "#abd9e9", "Negative" = "#E88B95")) +
  theme_void() + theme(legend.position = "none", strip.text = element_text(size = 12))
IcP_doughnut <- doughnut_data |> filter(Variable == "Ic_plank_slope") |> 
  ggplot(aes(x = 2, y = Fraction, fill = Sign)) + ggtitle("") +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + 
  scale_fill_manual(values = c("Positive" = "#abd9e9", "Negative" = "#E88B95")) +
  theme_void() + theme(legend.position = "none", strip.text = element_text(size = 12))
IcB_doughnut <- doughnut_data |> filter(Variable == "Ic_benthivorous_slope") |> 
  ggplot(aes(x = 2, y = Fraction, fill = Sign)) + ggtitle("") +
  geom_col(width = 1, color = "black") + coord_polar(theta = "y") + xlim(0.5, 2.5) + 
  scale_fill_manual(values = c("Positive" = "#abd9e9", "Negative" = "#E88B95")) +
  theme_void() + theme(legend.position = "none", strip.text = element_text(size = 12))

# Add pie plot insets
Figure_3A_inset <- Figure_3A +
  annotation_custom(grob = ggplotGrob(Gc_doughnut), xmin = 26, xmax = 41, ymin = 42, ymax = 48)
Figure_3B_inset <- Figure_3B+
  annotation_custom(grob = ggplotGrob(Fn_doughnut), xmin = 26, xmax = 41, ymin = 42, ymax = 48)
Figure_3C_inset <- Figure_3C +
  annotation_custom(grob = ggplotGrob(Fp_doughnut), xmin = 26, xmax = 41, ymin = 42, ymax = 48)
Figure_3D_inset <- Figure_3D +
  annotation_custom(grob = ggplotGrob(IcP_doughnut), xmin = 26, xmax = 41, ymin = 42, ymax = 48)
Figure_3E_inset <- Figure_3E +
  annotation_custom(grob = ggplotGrob(IcB_doughnut), xmin = 26, xmax = 41, ymin = 42, ymax = 48)

# Final Figure
Figure_3_tot = Figure_3A_inset / Figure_3B_inset / Figure_3C_inset / 
  Figure_3D_inset / Figure_3E_inset + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16))

#### Export the data  ----
## Figures
ggsave(Figure_3_tot, filename = "Figure_3.png", path = "Outputs/", device = "png", width = 6,  height = 14, dpi = 300)  