#### Setting up          ----
rm(list = ls())
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales") 
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") 
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("purrr")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("hexbin") ; library("ggforce")
library("RColorBrewer") ; library("MASS") ; library("sf") ; library("leaflet.extras")

## Functions
source("Scripts/00_functions_script.R")

## Download data
# Raw data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
TATB_WMED          <- read.delim("Data/TATB_WMED_1999-2021_clean.csv", sep = ";")
TATB_EMED          <- read.delim("Data/TATB_EMED_1999-2021_clean.csv", sep = ";")
Medits_total       <- rbind(TATB_WMED, TATB_EMED)
Nutrients_FishBase <- readxl::read_excel("Data/Nutrients_FishBase.xlsx") 
mdw_FishBase       <- readxl::read_excel("Data/mass_conversion_FishBase.xlsx") 
Metabolic_Rosen    <- read.delim("Data/Rosen_2025.csv", sep = ";")
Metabolic_Nina     <- read.delim("Data/Schiettekatte_2021_Metabolism.csv", sep = ",")
meow               <- st_read("Data/Datos Atlántico/MEOW-TNC/meow_ecos.shp")
load("Data/Datos Atlántico/community_and_traits.RData")
load("Data/Datos Atlántico/survey_all_combined.RData")
load("Data/Datos Atlántico/FishGlob_public_clean.RData")

### Exploration
###############

# Merging the three datasets
clean_numeric      <- function(x) as.numeric(gsub(",", ".", x))
data_Atl           <- data |> select(haul_id, longitude, latitude, year) |> distinct()
data_Med           <- Medits_total |> select(id, MEAN_LONGITUDE_DEC, MEAN_LATITUDE_DEC, YEAR) |> 
                      mutate(MEAN_LONGITUDE_DEC = clean_numeric(MEAN_LONGITUDE_DEC), 
                             MEAN_LATITUDE_DEC = clean_numeric(MEAN_LATITUDE_DEC)) |> distinct()
colnames(data_Med) <- colnames(data_Atl)
data_tot           <- rbind(data_Atl, data_Med) |> drop_na()

# Map projection
robin             <- "+proj=robin +datum=WGS84 +no_defs"
world             <- ne_countries(scale = "medium", returnclass = "sf")
world_robin       <- st_transform(world, crs = robin)
data_sf           <- st_as_sf(data_tot, coords = c("longitude", "latitude"), crs = 4326)
data_robin        <- st_transform(data_sf, crs = robin)
coords            <- st_coordinates(data_robin)
df_proj           <- data.frame(x = coords[,1], y = coords[,2])
df_proj           <- na.omit(df_proj)
bw                <- diff(range(df_proj$x)) / 250 # Cell size!
lon_seq           <- seq(-180, 180, by = 1)
lat_seq           <- seq(20, 90, by = 1)
coords            <- rbind(cbind(lon_seq, rep(20, length(lon_seq))), cbind(rep(180, length(lat_seq)), lat_seq), 
                           cbind(rev(lon_seq), rep(90, length(lon_seq))), cbind(rep(-180, length(lat_seq)), rev(lat_seq)))
coords            <- rbind(coords, coords[1, ])
poly_sf           <- st_sfc(st_polygon(list(coords)), crs = 4326)
poly_robin        <- st_transform(poly_sf, crs = robin)
world_union       <- st_transform(st_union(poly_robin), robin)
bb                <- st_as_sfc(st_bbox(world_union))
outside_robin     <- st_difference(bb, world_union)
world_robin_north <- world_robin[st_coordinates(st_centroid(world_robin))[ ,2] > 0, ]
lat_bounds        <- st_sfc(st_point(c(0, 20)), st_point(c(0, 90)), crs = 4326)
lat_bounds_robin  <- st_transform(lat_bounds, crs = robin)
y_lims            <- st_coordinates(lat_bounds_robin)[, 2]
coastline         <- world_robin |> dplyr::summarise()
lat20             <- st_sfc(st_linestring(cbind(seq(-180, 180, by = 1), rep(20, 361))), crs = 4326) |> st_transform(robin)
lat90             <- st_sfc(st_linestring(cbind(seq(-180, 180, by = 1), rep(90, 361))), crs = 4326) |> st_transform(robin)
lon180E           <- st_sfc(st_linestring(cbind(rep( 180, 111), seq(-20, 90, by = 1))), crs = 4326) |> st_transform(robin)
lon180W           <- st_sfc(st_linestring(cbind(rep(-180, 111), seq(-20, 90, by = 1))), crs = 4326) |> st_transform(robin)

# Long lat labels
lat_labels    <- data.frame(lon = c(-170, -170), lat = c(30, 40))
lat_labels_sf <- st_as_sf(lat_labels, coords = c("lon", "lat"), crs = 4326) |> st_transform(robin)
label_xy      <- st_coordinates(lat_labels_sf)
label_df      <- data.frame(x = label_xy[,1], y = label_xy[,2], label = c("30°N", "40°N"))

# Time series
data_sf    <- st_as_sf(data_tot, coords = c("longitude", "latitude"), crs = 4326)
data_robin <- st_transform(data_sf, crs = robin)
coords     <- st_coordinates(data_robin)
df_proj    <- data.frame(x = coords[,1], y = coords[,2])
hex_grid   <- st_make_grid(st_union(data_robin), cellsize = bw, square = FALSE)
hex_sf     <- st_sf(cell = seq_along(hex_grid), geometry = hex_grid)
pts_hex    <- st_join(data_robin, hex_sf, st_intersects)
hex_ts     <- pts_hex |> st_drop_geometry() |> group_by(cell, year) |> summarise(.groups = "drop") |>
  group_by(cell) |> summarise(ts = n()) |> mutate(ts_cat = case_when(ts < 5 ~ "<5", ts <= 10 ~ "5-10",
      ts <= 20 ~ "11-20", ts <= 30 ~ "21-30", ts <= 40 ~ "31-40", TRUE ~ ">40"),
      ts_cat = factor(ts_cat, levels = c("<5","5-10","11-20","21-30","31-40",">40"), ordered = TRUE))
hex_sf     <- hex_sf |> left_join(hex_ts, by = "cell")

# Regions
meow_robin <- st_transform(meow, robin)
hex_sf     <- st_join(hex_sf, meow_robin |> select(ECOREGION, ECO_CODE, PROVINCE, PROV_CODE))
hex_sf$ECOREGION[st_coordinates(st_centroid(hex_sf))[ ,1] > 13000000] <- "Aleutian Islands 2" 
hex_sf     <- hex_sf |> drop_na()
cell_meow  <- st_join(hex_sf, meow_robin)
eco_sf     <- hex_sf |> group_by(ECOREGION) |> summarise(geometry = st_union(geometry), .groups = "drop")
eco_pts    <- st_centroid(hex_sf)

# Plotting
(Figure_1a = ggplot() +
  geom_sf(data = poly_robin, fill = "lightblue", alpha = 0.4, color = NA) +
  geom_sf(data = hex_sf, aes(fill = ts_cat), color = "white", linewidth = 0.1) +
  geom_sf(data = coastline, fill = "grey90", color = "black", linewidth = 0.25) +
  geom_mark_rect(data = eco_pts, aes(x = st_coordinates(eco_pts)[,1], y = st_coordinates(eco_pts)[,2], group = ECOREGION), 
                 fill = "black", color = "black", alpha = .05, expand = unit(1.5, "mm")) +
  geom_sf(data = outside_robin, fill = "white", color = NA) +
  geom_sf(data = lat20, color = "black", linewidth = 1) +
  geom_sf(data = lat90, color = "black", linewidth = 1) +
  geom_sf(data = lon180E, color = "black", linewidth = 0.6) +
  geom_sf(data = lon180W, color = "black", linewidth = 0.6) +
  geom_text(data = label_df, aes(x = x, y = y, label = label), nudge_y = 3e5, size = 5.5, color = "grey30") +
  scale_fill_manual(name = "Number of years", values = c("<5" = "#1BDCF1", "5-10" = "#32A8E0", "11-20" = "#556CBC",
                                                         "21-30" = "#7343A9", "31-40" = "#9B0F8C", ">40" = "#DD1072"), drop = FALSE) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + 
  coord_sf(crs = st_crs(robin), ylim = y_lims, expand = FALSE) +
  theme_minimal() +
  theme(axis.title      = element_blank(),
        axis.text.y     = element_blank(),
        axis.text.x     = element_text(size = 16),
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "bottom"))

#####
## Merge properly both datasets
## MEDITS
MEDITS   <- rbind(TATB_WMED, TATB_EMED) |> data.frame() |> 
  dplyr::select(id, YEAR, MONTH, MEAN_LATITUDE_DEC, MEAN_LONGITUDE_DEC, HAUL_DURATION, SWEPT_AREA, MEAN_DEPTH, 
                BOTTOM_TEMPERATURE_BEGINNING, TOTAL_NUMBER_IN_THE_HAUL, TOTAL_WEIGHT_IN_THE_HAUL, MEDITS_CODE) |> 
  left_join(sp_code_list |> rename(MEDITS_CODE = spp.code)) |> 
  dplyr::filter(!str_detect(sci.name, "^NO\\s*")) |> dplyr::filter(catfau %in% c("Aa", "Ae", "Ao")) |> 
  rename(Species = sci.name) |> 
  mutate(Species = case_when(
    Species == "Trigloporus lastoviza" ~ "Chelidonichthys lastoviza",
    Species == "Pteromylaeus bovinus"  ~ "Aetomylaeus bovinus",
    Species == "Liza ramada"           ~ "Chelon ramada",
    Species == "Liza aurata"           ~ "Chelon auratus",
    Species == "Liza saliens"          ~ "Chelon saliens",
    TRUE ~ Species)) |> 
  left_join(rfishbase::load_taxa(), by = "Species") |> 
  mutate(Weight_CPUA = clean_numeric(TOTAL_WEIGHT_IN_THE_HAUL)/clean_numeric(SWEPT_AREA), 
         Num_CPUA = clean_numeric(TOTAL_NUMBER_IN_THE_HAUL)/clean_numeric(SWEPT_AREA)) |> 
  dplyr::select(-c(catfau, MEDITS_CODE, Subfamily, Order, Class, SuperClass, TOTAL_WEIGHT_IN_THE_HAUL, TOTAL_NUMBER_IN_THE_HAUL)) |> 
  mutate(Dataset = "MEDITS") |> 
  dplyr::select(Haul_ID = id, Year = YEAR, Month = MONTH, Latitude = MEAN_LATITUDE_DEC, Longitude = MEAN_LONGITUDE_DEC,
    Haul_duration = HAUL_DURATION, Haul_swept_area = SWEPT_AREA, Depth = MEAN_DEPTH, SBTemp = BOTTOM_TEMPERATURE_BEGINNING,
    SpecCode, Family, Genus, Species, Weight_CPUA, Num_CPUA, Dataset)

## FISHGLOB
FISHGLOB <- data |> data.frame() |> 
  mutate(haul_dur = case_when(is.na(haul_dur) ~ num / num_cpue, TRUE ~ haul_dur),
         haul_dur = case_when(is.na(haul_dur) ~ wgt / wgt_cpue, TRUE ~ haul_dur),
         area_swept = case_when(is.na(area_swept) ~ num / num_cpua, TRUE ~ area_swept),
         area_swept = case_when(is.na(area_swept) ~ wgt / wgt_cpua, TRUE ~ area_swept)) |> 
  dplyr::select(Haul_ID = haul_id, Year = year, Month = month, Latitude = latitude, Longitude = longitude, Haul_duration = haul_dur,
                Haul_swept_area = area_swept, Depth = depth, SBTemp = sbt, SpecCode, Family = family, Genus = genus, 
                Species = accepted_name, Weight_CPUA = wgt_cpua, Num_CPUA = num_cpua) |> 
  mutate(Dataset = "FISHGLOB", Haul_duration = if_else(Haul_duration > 1.62, Haul_duration / 60, Haul_duration),
         Haul_duration = Haul_duration*60, Weight_CPUA = Weight_CPUA*1000) # Work in grams

## FLUXGLOB
FLUXGLOB = rbind(MEDITS, FISHGLOB) |> 
  mutate(Longitude = as.numeric(str_replace(Longitude, ",", ".")),
         Latitude  = as.numeric(str_replace(Latitude, ",", "."))) |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |> st_transform(robin) |> 
  st_join(hex_sf |> select(cell, ECOREGION, ECO_CODE, PROVINCE, PROV_CODE), left = TRUE) |> drop_na(ECO_CODE)

# Code Figure 1 Ecoregions
Code_Fig_1a = data.frame(ECOREGION = unique(FLUXGLOB$ECOREGION),
                         ID_Fig1 = c(32, 33, 25, 34, 36, 35, 37, 38, 1, 3, 1, 4, 29, 27, 2, 24, 23, 10, 13, 12, 11, 15, 16, 
                                     5, 20, 21, 22, 18, 19, 17, 14, 30, 31, 28, 26, 7, 6, 8, 9))
Code_Fig_1b = data.frame(ID_Fig1 = seq(1,38,1),
                         subplot_Fig1 = c(rep("Fig1c", 9), rep("Fig1d", 13), rep("Fig1e", 4), 
                                          rep("Fig1f", 3), rep("Fig1g", 2), rep("Fig1h", 7)))
Code_Fig_1 = Code_Fig_1a |> left_join(Code_Fig_1b) 
FLUXGLOB = FLUXGLOB |> right_join(Code_Fig_1)

# Define the coverage area
Area_Covered = FLUXGLOB |> group_by(subplot_Fig1, cell) |> summarise(n = n()) |> group_by(subplot_Fig1) |> summarise(n = n()) |> 
  mutate(n = n * unique(round(as.numeric(st_area(hex_sf)) / 1e6, 2)))

### NMDS Communities
comm_eco_year         <- FLUXGLOB |> st_drop_geometry() |> group_by(ECOREGION, Year, SpecCode) |>
  summarise(Weight_CPUA = mean(Weight_CPUA, na.rm = TRUE),  .groups = "drop") |>
  pivot_wider(names_from = SpecCode, values_from = Weight_CPUA, values_fill = 0)
meta                  <- comm_eco_year |> select(ECOREGION, Year)
comm                  <- comm_eco_year |> select(-ECOREGION, -Year)
comm[is.na(comm)]     <- 0
comm_hel              <- vegan::decostand(comm, method = "hellinger")
keep                  <- rowSums(comm_hel) > 0
comm_hel              <- comm_hel[keep, ]
meta                  <- meta[keep, ]
nmds                  <- vegan::metaMDS(comm_hel, distance = "bray", k = 2, trymax = 100)
site_scores           <- as.data.frame(nmds$points)
colnames(site_scores) <- c("NMDS1", "NMDS2")
site_scores$ECOREGION <- meta$ECOREGION
site_scores$Year      <- meta$Year

### Kernel densities
k            <- MASS::kde2d(site_scores$NMDS1, site_scores$NMDS2, n = 200, lims = c(-4,3,-3,3.5))
z            <- as.vector(k$z)
dx           <- diff(k$x)[1]
dy           <- diff(k$y)[1]
p            <- z * dx * dy # probability mass per cell
ord          <- order(z, decreasing = TRUE)
cumprob      <- cumsum(p[ord]) / sum(p)
lev99        <- z[ord][which(cumprob >= 0.99)[1]]
lev75        <- z[ord][which(cumprob >= 0.75)[1]]
lev50        <- z[ord][which(cumprob >= 0.50)[1]]
data_isoline <- ggplot(site_scores, aes(NMDS1, NMDS2)) + stat_density_2d(geom = "contour", bins = 200) +
  scale_y_continuous(limits = c(-3,3.5)) + scale_x_continuous(limits = c(-4, 3))
isoline_0.99 <- ggplot_build(data_isoline)$data[[1]] |> dplyr::filter(level == 0.00135)
isoline_0.75 <- ggplot_build(data_isoline)$data[[1]] |> dplyr::filter(level == 0.05805)
isoline_0.50 <- ggplot_build(data_isoline)$data[[1]] |> dplyr::filter(level == 0.10665)

# Plotting
(Figure_1b <- ggplot(site_scores, aes(NMDS1, NMDS2)) +
    geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = "#ECFFDC", color = "black", linewidth = 1) +
    geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = "#C1E1C1", color = "black") +
    geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = "#93C572", color = "black") +
    geom_point(fill = "#2E8B57", color = "#B4C424", size = 1.5, shape = 21, alpha = .7) +
    scale_y_continuous(name = "NMDS2", limits = c(-3,3.5)) + scale_x_continuous(name = "NMDS1", limits = c(-4, 3)) +
    theme_classic() +
    theme(axis.title      = element_text(size = 18),
          axis.text       = element_text(size = 16),
          legend.position = "none",
          panel.border    = element_rect(color = "black", fill = NA, linewidth = 1)))

# Explore dominant species for ER and Year
sites_Kernel <- site_scores |> dplyr::filter(NMDS1 > 0.19, NMDS1 < 0.53, NMDS2 > 0.22, NMDS2 < 1.01) |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = ""))
sites_Kernel <- comm_eco_year |> mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% sites_Kernel$ER_Year) |> 
  pivot_longer(cols = -c(ECOREGION, Year, ER_Year), names_to = "SpecCode", values_to = "Biomass") |> 
  mutate(SpecCode = as.integer(SpecCode)) |> left_join(rfishbase::load_taxa() |> select(SpecCode, Species)) |> 
  group_by(SpecCode, Species) |> summarise(Biomass = mean(Biomass))

## Sub-NMDS Fig1c-1h
# Color gradient
year_breaks <- c(1963, seq(1970, 2020, 10), 2024)
year_cols <- c("#FCEA90","#F0D90E","#F0BF0E","#F0A50E","#E38C09","#E36709","#E33109","#C42D0A")

# Figure 1c
Fig1c_sites = FLUXGLOB |> dplyr::filter(subplot_Fig1 == "Fig1c") |> 
  mutate(ECOREGION = case_match(ECOREGION, "Aleutian Islands 2" ~ "Aleutian Islands", .default = ECOREGION)) |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> distinct(ER_Year)
site_scores_1c = site_scores |> 
  mutate(ECOREGION = case_match(ECOREGION, "Aleutian Islands 2" ~ "Aleutian Islands", .default = ECOREGION)) |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% Fig1c_sites$ER_Year) |> dplyr::select(-ER_Year)

Figure_1c = ggplot(site_scores_1c, aes(NMDS1, NMDS2)) +
  geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = NA, color = "black", linewidth = 0.75) +
  geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_point(aes(fill = Year, color = Year), size = 2.5, shape = 21) +
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_color_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_y_continuous(name = "", limits = c(-2.5,3.5)) + scale_x_continuous(name = "", limits = c(-4, 3)) +
  theme_void() +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",
        panel.border    = element_blank())

# Figure 1d
Fig1d_sites = FLUXGLOB |> dplyr::filter(subplot_Fig1 == "Fig1d") |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> distinct(ER_Year)
site_scores_1d = site_scores |> mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% Fig1d_sites$ER_Year) |> dplyr::select(-ER_Year)

Figure_1d = ggplot(site_scores_1d, aes(NMDS1, NMDS2)) +
  geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = NA, color = "black", linewidth = 0.75) +
  geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_point(aes(fill = Year, color = Year), size = 2.5, shape = 21) +
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_color_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_y_continuous(name = "", limits = c(-2.5,3.5)) + scale_x_continuous(name = "", limits = c(-4, 3)) +
  theme_void() +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",
        panel.border    = element_blank())

# Figure 1e
Fig1e_sites = FLUXGLOB |> dplyr::filter(subplot_Fig1 == "Fig1e") |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> distinct(ER_Year)
site_scores_1e = site_scores |> mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% Fig1e_sites$ER_Year) |> dplyr::select(-ER_Year)

Figure_1e = ggplot(site_scores_1e, aes(NMDS1, NMDS2)) +
  geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = NA, color = "black", linewidth = 0.75) +
  geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_point(aes(fill = Year, color = Year), size = 2.5, shape = 21) +
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_color_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_y_continuous(name = "", limits = c(-2.5,3.5)) + scale_x_continuous(name = "", limits = c(-4, 3)) +
  theme_void() +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",
        panel.border    = element_blank())

# Figure 1f
Fig1f_sites = FLUXGLOB |> dplyr::filter(subplot_Fig1 == "Fig1f") |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> distinct(ER_Year)
site_scores_1f = site_scores |> mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% Fig1f_sites$ER_Year) |> dplyr::select(-ER_Year)

Figure_1f = ggplot(site_scores_1f, aes(NMDS1, NMDS2)) +
  geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = NA, color = "black", linewidth = 0.75) +
  geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_point(aes(fill = Year, color = Year), size = 2.5, shape = 21) +
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_color_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_y_continuous(name = "", limits = c(-2.5,3.5)) + scale_x_continuous(name = "", limits = c(-4, 3)) +
  theme_void() +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",
        panel.border    = element_blank())

# Figure 1g
Fig1g_sites = FLUXGLOB |> dplyr::filter(subplot_Fig1 == "Fig1g") |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> distinct(ER_Year)
site_scores_1g = site_scores |> mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% Fig1g_sites$ER_Year) |> dplyr::select(-ER_Year)

Figure_1g = ggplot(site_scores_1g, aes(NMDS1, NMDS2)) +
  geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = NA, color = "black", linewidth = 0.75) +
  geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_point(aes(fill = Year, color = Year), size = 2.5, shape = 21) +
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_color_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_y_continuous(name = "", limits = c(-2.5,3.5)) + scale_x_continuous(name = "", limits = c(-4, 3)) +
  theme_void() +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",
        panel.border    = element_blank())

# Figure 1h
Fig1h_sites = FLUXGLOB |> dplyr::filter(subplot_Fig1 == "Fig1h") |> 
  mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> distinct(ER_Year)
site_scores_1h = site_scores |> mutate(ER_Year = paste(ECOREGION, "_", Year, sep = "")) |> 
  dplyr::filter(ER_Year %in% Fig1h_sites$ER_Year) |> dplyr::select(-ER_Year)

Figure_1h = ggplot(site_scores_1h, aes(NMDS1, NMDS2)) +
  geom_polygon(data = isoline_0.99, aes(x, y, group = group), fill = NA, color = "black", linewidth = 0.75) +
  geom_polygon(data = isoline_0.75, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_polygon(data = isoline_0.50, aes(x, y, group = group), fill = NA, color = "black", linetype = "dashed") +
  geom_point(aes(fill = Year, color = Year), size = 2.5, shape = 21) +  
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_color_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024), breaks = year_breaks) +
  scale_y_continuous(name = "", limits = c(-2.5,3.5)) + scale_x_continuous(name = "", limits = c(-4, 3)) +
  theme_void() +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        legend.position = "none",
        panel.border    = element_blank())

## Color gradient legend
heatmap_df <- tibble(Year = 1963:2024, y = 1)

Figure_1i <- ggplot(heatmap_df, aes(Year, y, fill = Year)) + geom_tile(show.legend = F) + 
  scale_fill_gradientn(colours = year_cols, values = scales::rescale(year_breaks), limits = c(1963, 2024)) +
  scale_x_continuous(breaks = year_breaks, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  theme_void() +
  theme(axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        panel.grid   = element_blank(),
        axis.title   = element_blank(),
        panel.border = element_rect(color = "white", fill = NA, linewidth = 1))

##### Export Figure
ggsave(Figure_1a, filename = "Outputs/FLUXGLOB/Raw/Figure_1a_Eli.png", path = "Outputs/", device = "png", width = 16, height = 04, dpi = 300)  
ggsave(Figure_1b, filename = "Outputs/FLUXGLOB/Raw/Figure_1b_Eli.png", path = "Outputs/", device = "png", width = 07, height = 07, dpi = 300)  
ggsave(Figure_1c, filename = "Outputs/FLUXGLOB/Raw/Figure_1c_Eli.png", path = "Outputs/", device = "png", width = 02, height = 02, dpi = 300)  
ggsave(Figure_1d, filename = "Outputs/FLUXGLOB/Raw/Figure_1d_Eli.png", path = "Outputs/", device = "png", width = 02, height = 02, dpi = 300)  
ggsave(Figure_1e, filename = "Outputs/FLUXGLOB/Raw/Figure_1e_Eli.png", path = "Outputs/", device = "png", width = 02, height = 02, dpi = 300)  
ggsave(Figure_1f, filename = "Outputs/FLUXGLOB/Raw/Figure_1f_Eli.png", path = "Outputs/", device = "png", width = 02, height = 02, dpi = 300)  
ggsave(Figure_1g, filename = "Outputs/FLUXGLOB/Raw/Figure_1g_Eli.png", path = "Outputs/", device = "png", width = 02, height = 02, dpi = 300)  
ggsave(Figure_1h, filename = "Outputs/FLUXGLOB/Raw/Figure_1h_Eli.png", path = "Outputs/", device = "png", width = 02, height = 02, dpi = 300)  
ggsave(Figure_1i, filename = "Outputs/FLUXGLOB/Raw/Figure_1i_Eli.png", path = "Outputs/", device = "png", width = 10, height = .2, dpi = 300)  