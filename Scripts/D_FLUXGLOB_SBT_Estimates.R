#### Setting up     ----
rm(list = ls())
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales") ; library("data.table")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("ggforce") ; library("ggtext")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("purrr") ; library("sf")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("hexbin") ; library("leaflet.extras") ; library("brms")
library("RColorBrewer") ; library("MASS") ; library("stringr") ; library("terra") ; library("reticulate")

## Data & Functions
source("Scripts/00_functions_script.R")
source("Scripts/A_FLUXGLOB_Context.R")
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB.RData")

# ---- DATA DOWNLOAD ----
reticulate::py_config() ; reticulate::py_require(c("copernicusmarine")) ; cm <- import("copernicusmarine")
cm$login(username = Sys.getenv("COPERNICUSMARINE_USERNAME"), password = Sys.getenv("COPERNICUSMARINE_PASSWORD"))
cm$subset(dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1M-m",
  variables          = list("bottomT"),
  minimum_longitude  = -180,
  maximum_longitude  = 180,
  minimum_latitude   = 20,
  maximum_latitude   = 90,
  start_datetime     = "1993-01-01T00:00:00",
  end_datetime       = "2025-12-31T00:00:00",
  output_directory   = "data",
  output_filename    = "bottomT_1993_2025.nc")
bt      <- terra::rast("data/bottomT_1993_2025.nc")
months  <- as.integer(format(time(bt), "%m"))
bt_clim <- tapp(bt, index = months, fun = mean, na.rm = TRUE)

# One record per haul
hauls         <- FLUXGLOB_cor |> st_as_sf() |> distinct(Haul_ID, Month, .keep_all = TRUE) |> st_transform(4326)
hauls_vect    <- vect(hauls)
hauls$Month   <- as.integer(hauls$Month)
hauls$botTemp <- NA_real_
for (m in 1:12) {idx <- which(hauls$Month == m)
  if (length(idx) > 0) { hauls$botTemp[idx] <- terra::extract(bt_clim[[m]], hauls_vect[idx])[, 2]}}
FLUXGLOB_cor  <- FLUXGLOB_cor |> left_join(st_drop_geometry(hauls[, c("Haul_ID", "botTemp")]), by = "Haul_ID")

## Visualization
# Convert climatology raster to points
haul_cells   <- FLUXGLOB_cor |> distinct(cell, Month) |>
  mutate(Month = factor(month.abb[as.integer(Month)], levels = month.abb, ordered = TRUE))
bt_df        <- as.data.frame(bt_clim, xy = TRUE, na.rm = TRUE)
names(bt_df) <- c("longitude", "latitude", month.abb)
bt_long      <- bt_df |> pivot_longer(Jan:Dec, names_to = "Month", values_to = "bottomT") |>
  mutate(Month = factor(Month, levels = month.abb, ordered = TRUE))
bt_sf        <- st_as_sf(bt_long, coords = c("longitude", "latitude"), crs = 4326)
bt_robin     <- st_transform(bt_sf, robin)
bt_hex       <- st_join(bt_robin, hex_sf["cell"])
hex_temp     <- bt_hex |> st_drop_geometry() |> group_by(cell, Month) |> 
  summarise(bottomT = mean(bottomT, na.rm = T), .groups = "drop") |>
  inner_join(haul_cells, by = c("cell", "Month"))
hex_plot     <- left_join(hex_sf, hex_temp, by = "cell", relationship = "many-to-many") |> 
  inner_join(haul_cells, by = c("cell", "Month"), relationship = "many-to-many") |> drop_na(Month)

(Figure_S2 = ggplot() + geom_sf(data = poly_robin, fill = "lightblue", colour = NA) +
  geom_sf(data = hex_plot, aes(fill = bottomT, colour = bottomT), linewidth = 0.1) +
  geom_sf(data = coastline, fill = "grey90", colour = "black", linewidth = 0.25) +
  geom_sf(data = outside_robin, fill = "white", color = NA) +
  geom_sf(data = lat20, color = "black", linewidth = 1) +
  geom_sf(data = lat90, color = "black", linewidth = 1) +
  geom_sf(data = lon180E, color = "black", linewidth = 0.6) +
  geom_sf(data = lon180W, color = "black", linewidth = 0.6) +
  scale_fill_viridis_c(name = "Sea Bottom temperature (°C)", na.value = "white", option = "plasma") +
  scale_color_viridis_c(name = "Sea Bottom temperature (°C)", na.value = "white", option = "plasma") +
  coord_sf(crs = st_crs(robin), ylim = y_lims, expand = FALSE ) +
  theme_minimal() + facet_wrap(~Month, ncol = 2) +
  theme(plot.margin = margin(0, 0, 0, 0),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom"))

#### Export the data  ----
save(FLUXGLOB_cor, file = "Outputs/FLUXGLOB/dat_proc/FLUXGLOB_with_temp.RData")

##### Export Figure   ----
ggsave(Figure_S2, filename = "FLUXGLOB/Raw/Figure_S2_Eli.png", path = "Outputs/", device = "png", width = 10, height = 8, dpi = 300)