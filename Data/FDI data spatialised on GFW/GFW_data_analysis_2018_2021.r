rm(list = ls(all = TRUE)) # Clear the entire workspace, removing all objects

library(dplyr)            # Load dplyr for data manipulation
library(ggplot2)          # Load ggplot2 for plotting
library(sf)               # Load sf for handling spatial vector data
library(rnaturalearth)    # Load rnaturalearth for country boundaries
library(rnaturalearthdata)# Additional natural earth data
library(raster)           # For handling raster data

wd <- "D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/B-USEFUL/Task 4.2/MED_ANALYSIS/data/GFW" # Set working directory path
setwd(wd)  # Change R's current working directory
resdir <- "D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/B-USEFUL/Task 4.2/MED_ANALYSIS/data/GFW/results_gfw" # Directory for results

# Define GSA codes
gsas <- c("1","2","5","6","7","8","9","10","11","15","16","17","18","19",
          "20","21","22","23","25","11.1","11,2")
gsas <- paste("GSA", gsas, sep="") # Prepend "GSA" to each code

fleet <- "demersal fisheries" # Could be "demersal fisheries", "demersal+pelagic", or "trawlers"
countries_codes <- c("ESP", "FRA", "ITA", "MLT", "HRV", "SVN", "GRC","CYP") # Codes for countries
countries <- c("Spain", "France", "Italy", "Malta", "Croatia", "Slovenia", "Greece","Cyprus") # Country names
tab_country <- data.frame(country_code = countries_codes, country = countries) # Create lookup table

fdi_years <- ys <- c(2018:2021) # Define the range of years for analysis
source("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/SEAwise/5.3/degrees_to_csquare.r") # Convert lat-lon to csquare
source("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/SEAwise/5.3/csquare_2_latlon.r")   # Convert csquare back to lat-lon
depth <- raster::raster("D:/GIS/EMODNET_bathymetry/EMODNET_bathymetry_2022.tif") # Raster of bathymetry data

files <- list() # Initialize a list to store fishing effort data
files[[1]] <- read.table("2021_public-global-fishing-effort-v3.0.csv", sep=",", head=TRUE) # 2021 data
files[[1]]$year <- 2021
files[[2]] <- read.table("2022_public-global-fishing-effort-v3.0.csv", sep=",", head=TRUE) # 2022 data
files[[2]]$year <- 2022
files[[3]] <- read.table("2023_public-global-fishing-effort-v3.0.csv", sep=",", head=TRUE) # 2023 data
files[[3]]$year <- 2023
# The following lines are commented out but would handle 2024 data if needed:
# files[[4]] <- read.table("2024_public-global-fishing-effort-v3.0.csv",sep=",", head=TRUE)
# files[[4]]$year <- 2024

d <- do.call(rbind,files)        # Combine all data frames into one
d <- d[d$Flag %in% countries_codes, ] # Filter rows to keep only selected countries

#-----------
fdi <- read.table("FDI Effort by country.csv", sep=";", header=TRUE) # Load FDI effort data
fdi <- fdi[fdi$Country %in% countries & fdi$Year %in% fdi_years & fdi$Confidential =="N", ] # Filter for chosen countries, years, and non-confidential data
fdi <- fdi[fdi$Vessel.Length.Category %in% c("VL1218","VL1824","VL2440","VL40XX"), ]        # Filter by vessel length categories
fdi <- fdi[fdi$Sub.region %in% gsas, ]  # Filter for specified GSAs
#-----------

# Determine the gear list based on the chosen fleet
if (fleet == "demersal fisheries") {
  
  # Define a list of gears for demersal fisheries
  gear_list <- list(
    dredge_fishing = c("DRB"),
    fixed_gear = c("GTR"),
    pots_and_traps = c("FPO", "FYK"),
    seiners = c("SB","SV","SPR"),
    set_gillnets = c("GNS"),
    set_longlines = c("LLS"),
    squid_jigger = c("LHP", "LHM"),
    trawlers = c("OTB", "OTT", "TBB", "PTB")
  )
  
  # Create a long-format data frame from the list
  gear_table <- do.call(rbind, lapply(names(gear_list), function(gear) {
    data.frame(gear = gear, contents = gear_list[[gear]], stringsAsFactors = FALSE)
  }))
  
} else if (fleet == "demersal+pelagic"){
  
  # Two gear lists for demersal and pelagic
  gear_list1 <- list(
    dredge_fishing = c("DRB"),
    fixed_gear = c("GTR"),
    pots_and_traps = c("FPO", "FYK"),
    seiners = c("SB", "SV", "SPR"),
    set_gillnets = c("GNS"),
    set_longlines = c("LLS"),
    squid_jigger = c("LHP", "LHM"),
    trawlers = c("OTB", "OTT", "TBB", "PTB")
  )
  
  gear_list2 <- list(
    drifting_longlines = c("LLD"),
    fixed_gear = c("GND"),
    pots_and_traps = c("FPN"),
    purse_seines = c("PS", "LA"),
    trawlers = c("OTM", "PTM"),
    trollers = c("LTL")
  )
  
  # Merge both lists and remove duplicates
  merged_gear_list <- lapply(
    split(unlist(c(gear_list1, gear_list2)),
          rep(names(c(gear_list1, gear_list2)), lengths(c(gear_list1, gear_list2)))),
    unique
  )
  
  # Create a long-format data frame from merged gear lists
  gear_table <- do.call(rbind, lapply(names(merged_gear_list), function(gear) {
    data.frame(gear = gear, contents = merged_gear_list[[gear]], stringsAsFactors = FALSE)
  }))
  
  print(gear_table) # Check the gear table
  
} else if (fleet == "trawlers"){
  
  # Define a list of gears for trawlers
  gear_list <- list(
    trawlers = c("OTB", "OTT", "TBB", "PTB")
  )
  
  # Create a long-format data frame
  gear_table <- do.call(rbind, lapply(names(gear_list), function(gear) {
    data.frame(gear = gear, contents = gear_list[[gear]], stringsAsFactors = FALSE)
  }))
}

#------------
# Adjust the coordinates by a small offset
d$Lon <- d$Lon + 0.05
d$Lat <- d$Lat + 0.05

# Convert lat/lon to csquares
d$csquares <- degrees_to_csquare(
  data = d,
  grid_square = 0.1,
  latitude_name = "Lat",
  longitude_name = "Lon",
  boundary_ajustement_factor = 1e-06
)[, ncol(d)+1]

# Retrieve the lat/lon from csquares
coord <- csquare_2_latlon(d$csquares)
d$Lat <- coord[[1]]
d$Lon <- coord[[2]]
#------------

# Filter FDI data to keep only the gears in gear_table
fdi <- fdi[fdi$Gear.Type %in% gear_table$contents, ]

# Merge FDI with gear_table to unify gear type info
fdi2 <- fdi %>%
  left_join(gear_table, by = c("Gear.Type" = "contents"))
fdi <- data.frame(fdi2)
fdi$Total.Fishing.Days <- as.numeric(fdi$Total.Fishing.Days) # Ensure effort is numeric

# Summarize fishing days by year, country, and gear
tab <- fdi %>%
  dplyr::group_by(Year, Country, gear) %>%
  dplyr::summarise(effort = sum(Total.Fishing.Days, na.rm=TRUE))
tab$code <- paste(tab$Country, tab$gear, sep="_") # Create an identifier

#### ---- Effort spatial pattern ----
# Filter GFW data to keep only the relevant gear
d <- d[d$Geartype %in% unique(gear_table$gear), ]

# Group by multiple columns and sum the Apparent.Fishing.Hours
fp <- d %>%
  dplyr::group_by(csquares, year, Lat, Lon, Flag, Geartype) %>%
  dplyr::summarise(effort = sum(Apparent.Fishing.Hours, na.rm=TRUE))

# For each csquare, Flag, and Geartype, take the mean effort
fp <- fp %>%
  dplyr::group_by(csquares, Lat, Lon, Flag, Geartype) %>%
  dplyr::summarise(effort = mean(effort, na.rm=TRUE))

# Merge with country lookup to replace Flag codes
fp <- fp %>%
  left_join(tab_country, by = c("Flag" = "country_code"))
fp$code <- paste(fp$country, fp$Geartype, sep="_") # Create an identifier
fp <- data.frame(fp)

# Calculate the proportion of effort for each code
fp_st <- fp %>% 
  dplyr::group_by(code) %>%
  dplyr::mutate(prop = effort / sum(effort, na.rm=TRUE)) %>%
  dplyr::ungroup()

fp <- as.data.frame(fp_st)

#-----------
# Identify codes that appear in both FDI data and GFW data
fp_codes <- unique(fp$code)
fdi_codes <- unique(tab$code)
fp_codes
fdi_codes
fdi_codes[which(!fdi_codes %in% fp_codes )]
fp_codes[which(!fp_codes %in% fdi_codes )]

common_codes <- intersect(fp_codes, fdi_codes)
print(common_codes)

# Filter both data sources to keep only the intersection
tab <- tab[tab$code %in% common_codes, ]
fp <- fp[fp$code %in% common_codes, ]

# #-----------
# # Below are lines for creating a map with facet_grid, commented out for reference
# world <- ne_countries(scale = "medium", returnclass = "sf")
# med_extent <- c(xmin = -6, xmax = 37, ymin = 30, ymax = 48)
# plot_effort <- ggplot() +
#     geom_sf(data = world, fill = "gray80", color = "black") +
#     geom_tile(data = fp[!is.na(fp$effort),], aes(x = Lon, y = Lat, fill = effort)) +
#     scale_fill_viridis_c(option = "plasma", na.value = "white") +
#     facet_grid(rows = vars(Flag), cols = vars(Geartype)) +
#     coord_sf(xlim = c(med_extent["xmin"], med_extent["xmax"]), 
#              ylim = c(med_extent["ymin"], med_extent["ymax"]), expand = FALSE) +
#     theme_minimal() +
#     labs(title = "Fishing Effort by Country and Gear Type",
#          x = "Longitude",
#          y = "Latitude",
#          fill = "Effort")
# print(plot_effort)
# #---------------

y = 1          # Initialize a counter (used in the loop)
exp <- list()  # Create an empty list to store temporary data frames

# Loop over the FDI years
for (y in 1:length(ys)) {
  
  temp <- fp
  colnames(temp)[which(colnames(temp)=="effort")] <- "effort_gfw" # Rename effort column to avoid confusion
  taby <- tab[tab$Year == ys[y], c("Year","effort","code")]       # Filter tab by current year
  temp <- temp %>% left_join(taby, by = c("code" = "code"))       # Merge with GFW data
  temp <- temp[!is.na(temp$effort), ]                             # Remove rows with no matching FDI data
  temp$FD <- temp$prop * temp$effort                              # Multiply total effort by proportion
  exp[[length(exp)+1]] <- temp                                    # Store results in the list
  rm(temp)                                                        # Remove temporary variable
}

# Combine all yearly results
dd <- do.call(rbind, exp)
dd <- data.frame(dd)
head(dd) # Check first rows

# Group by csquare, lat, lon, year, and sum FD
dd_aggregated <- dd %>%
  dplyr::group_by(csquares, Lat, Lon, Year) %>%
  dplyr::summarise(FD_total = sum(FD, na.rm=TRUE)) %>%
  dplyr::ungroup()

print(head(dd_aggregated)) # Inspect the aggregated data

# Extract depth from the raster using coordinates
dd_aggregated$depth <- raster::extract(depth, dd_aggregated[, c("Lon","Lat")])
dd_aggregated$depth  <- -1 * dd_aggregated$depth  # Convert depth to positive (if negative)
dd_aggregated <- dd_aggregated[dd_aggregated$depth <= 1000, ] # Keep only depths <= 1000m
dd_aggregated <- dd_aggregated[!is.na(dd_aggregated$Year) & dd_aggregated$FD_total > 0,] 
# Remove records with no year or zero effort

# Summarize total FD per year
dd_aggregated2 <- dd_aggregated %>%
  group_by(Year) %>%
  summarise(FD = sum(FD_total, na.rm=TRUE)) %>%
  filter(!is.na(Year))

# Plot regression of FD over years
plot_regression <- ggplot(dd_aggregated2, aes(x = Year, y = FD)) +
  geom_point() +  
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = TRUE) +
  theme_minimal() +
  labs(title = "Linear Regression of FD over Years",
       x = "Year",
       y = "FD (Total)")

print(plot_regression) # Display the regression plot

# Spearman correlation
cor.test(dd_aggregated2$Year, dd_aggregated2$FD, method = "spearman")

# The following lines (commented out) would add 2012 and 2013 data if needed:
# d2013 <- dd_aggregated[dd_aggregated$Year == 2013, ]
# d2012 <- d2013
# d2012$Year <- 2012
# dd_aggregated <- rbind(dd_aggregated, d2012)

dd_aggregated$fleet <- fleet # Tag with the fleet name

# Load world map (country boundaries)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the bounding box for the Mediterranean Sea
med_extent <- c(xmin = -6, xmax = 37, ymin = 30, ymax = 48)

# Prepare label for the map
label <- paste("Fishing Effort of", fleet)

# Create a facet_wrap map by year
plot_effort <- ggplot() +
  theme_minimal() +
  geom_tile(data = dd_aggregated, aes(x = Lon, y = Lat, fill = FD_total)) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = "white", trans="pseudo_log") +
  facet_wrap(.~Year) +
  geom_sf(data = world, fill = "gray80", color = "black") +
  coord_sf(xlim = c(med_extent["xmin"], med_extent["xmax"]),
           ylim = c(med_extent["ymin"], med_extent["ymax"]), expand = FALSE) +
  labs(title = label,
       x = "Longitude",
       y = "Latitude",
       fill = "Effort")

# Build the output file name for saving the plot
output_file <- file.path(resdir, paste0("GFW_FDI_rescaled_effort_", paste0(range(ys), collapse="_"), ".jpg"))

# Save the plot to a file
ggsave(filename = output_file, plot = plot_effort, width = 10, height = 6, dpi = 300)

#---------------

# Write the aggregated data to a CSV file
write.table(dd_aggregated, file.path(resdir, paste0("GFW_FDI_rescaled_effort_",paste0(range(ys),collapse="_"),".csv")),
            sep=";", row.names=FALSE)

# Create a subfolder named 'rasters' to store raster files
dir.create(file.path(resdir,"rasters"))

# Create and export raster files for each year
for (year in unique(dd_aggregated$Year)) {
  
  # Filter the data for the current year
  dd_year <- dd_aggregated %>% filter(Year == year)
  
  # Convert the data to a raster
  r <- rasterFromXYZ(dd_year[, c("Lon", "Lat", "FD_total")])
  
  # Set the CRS to WGS84
  crs(r) <- CRS("+proj=longlat +datum=WGS84")
  
  # Construct the raster file name
  output_file <- file.path(resdir,"rasters", paste0(year,"_FE_",fleet, ".tif"))
  
  # Write the raster to disk
  writeRaster(r, output_file, format = "GTiff", overwrite = TRUE)
  
  message(paste("Raster saved for year:", year)) # Inform about the saved raster
}

