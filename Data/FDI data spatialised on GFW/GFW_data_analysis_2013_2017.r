rm(list = ls(all = TRUE)) # Clear the entire workspace, removing all objects

library(dplyr)            # Load the dplyr package for data manipulation
library(ggplot2)          # Load the ggplot2 package for data visualization
library(sf)               # Load sf for handling geospatial vector data
library(rnaturalearth)    # Load rnaturalearth for world map data
library(rnaturalearthdata)# Additional natural earth map data
library(raster)           # Handle raster (grid) data

wd <- "D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/B-USEFUL/Task 4.2/MED_ANALYSIS/data/GFW" # Set working directory path
setwd(wd)  # Set the working directory
resdir <- "D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/B-USEFUL/Task 4.2/MED_ANALYSIS/data/GFW/results_gfw" # Define results directory

gsas <- c("1","2","5","6","7","8","9","10","11","15","16","17","18","19",
          "20","21","22","23","25","11.1","11,2")  # Vector of GSA codes
gsas <- paste("GSA",gsas,sep="") # Combine "GSA" prefix with each code

fleet <- "demersal fisheries" # Define the fleet type; alternatives: "demersal fisheries", "demersal+pelagic", "trawlers"
countries_codes <- c("ESP", "FRA", "ITA", "MLT", "HRV", "SVN", "GRC","CYP") # Codes for countries
countries <- c("Spain", "France", "Italy", "Malta", "Croatia", "Slovenia", "Greece","Cyprus") # Country names
tab_country <- data.frame(country_code = countries_codes, country = countries) # Create a lookup table for codes

fdi_years <- ys <- c(2013:2017) # Define the range of years to analyze
source("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/SEAwise/5.3/degrees_to_csquare.r") # Custom function for csquare conversion
source("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/SEAwise/5.3/csquare_2_latlon.r")   # Custom function to get lat/lon from csquare
depth <- raster::raster("D:/GIS/EMODNET_bathymetry/EMODNET_bathymetry_2022.tif") # Raster layer for bathymetry data

files <- list() # Initialize a list to store data files
files[[1]] <- read.table("2016_public-global-fishing-effort-v3.0.csv",sep=",", head=TRUE) # Read the 2016 data
files[[1]]$year <- 2016 # Assign year 2016
files[[2]] <- read.table("2017_public-global-fishing-effort-v3.0.csv",sep=",", head=TRUE) # Read the 2017 data
files[[2]]$year <- 2017 # Assign year 2017
# files[[3]] <- read.table("2023_public-global-fishing-effort-v3.0.csv",sep=",", head=TRUE)
# files[[3]]$year <- 2023
# files[[4]] <- read.table("2024_public-global-fishing-effort-v3.0.csv",sep=",", head=TRUE)
# files[[4]]$year <- 2024
# The above lines are commented out; they would read data for 2023 and 2024 if needed

d <- do.call(rbind,files)    # Combine the data frames in 'files' into one
d <- d[d$Flag %in% countries_codes, ] # Filter rows where Flag matches the selected countries

#-----------
fdi <- read.table("FDI Effort by country.csv", sep=";",header=TRUE) # Read FDI effort data
fdi <- fdi[fdi$Country %in% countries & fdi$Year %in% fdi_years & fdi$Confidential =="N", ] # Filter on countries, years, and non-confidential
fdi <- fdi[fdi$Vessel.Length.Category %in% c("VL1218","VL1824","VL2440","VL40XX"), ]         # Filter on specific vessel length categories
fdi <- fdi[fdi$Sub.region %in% gsas, ] # Filter on the specified GSA regions
#-----------

if (fleet == "demersal fisheries") {
  
  # Definition of gear list for demersal fisheries
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
  
  # Create a long-format data frame from gear_list
  gear_table <- do.call(rbind, lapply(names(gear_list), function(gear) {
    data.frame(gear = gear, contents = gear_list[[gear]], stringsAsFactors = FALSE)
  }))
  
}  else if (fleet == "demersal+pelagic"){
  
  # Two gear lists, one for demersal and one for pelagic
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
  
  # Merge the two lists, removing duplicates
  merged_gear_list <- lapply(
    split(unlist(c(gear_list1, gear_list2)), 
          rep(names(c(gear_list1, gear_list2)), lengths(c(gear_list1, gear_list2)))),
    unique
  )
  
  # Create a long-format data frame from the merged list
  gear_table <- do.call(rbind, lapply(names(merged_gear_list), function(gear) {
    data.frame(gear = gear, contents = merged_gear_list[[gear]], stringsAsFactors = FALSE)
  }))
  
  print(gear_table) # Print the table for verification
  
} else if (fleet == "trawlers"){
  
  # Definition of gear list for trawlers only
  gear_list <- list(
    trawlers = c("OTB", "OTT", "TBB", "PTB")
  )
  
  # Create a long-format data frame from gear_list
  gear_table <- do.call(rbind, lapply(names(gear_list), function(gear) {
    data.frame(gear = gear, contents = gear_list[[gear]], stringsAsFactors = FALSE)
  }))
}

#------------
d$Lon <- d$Lon + 0.05  # Adjust longitude by a small offset
d$Lat <- d$Lat + 0.05  # Adjust latitude by a small offset
d$csquares <- degrees_to_csquare(data=d, grid_square=0.1, latitude_name="Lat", longitude_name="Lon", boundary_ajustement_factor = 1e-06)[,ncol(d)+1] 
# Convert lat/lon to "csquares" identifier

coord <- csquare_2_latlon(d$csquares) # Retrieve lat/lon from csquare
d$Lat <- coord[[1]] # Overwrite the latitude with csquare-derived lat
d$Lon <- coord[[2]] # Overwrite the longitude with csquare-derived lon
#------------

fdi <- fdi[fdi$Gear.Type %in% gear_table$contents, ] # Filter the FDI data to include only the gears defined in gear_table

fdi2 <- fdi %>%
  left_join(gear_table, by = c("Gear.Type" = "contents")) # Merge to associate gear type categories
fdi <- data.frame(fdi2) # Convert to a standard data frame
fdi$Total.Fishing.Days <- as.numeric(fdi$Total.Fishing.Days) # Ensure effort is numeric

tab <- fdi %>% dplyr::group_by(Year,Country,gear) %>% dplyr::summarise(effort = sum(Total.Fishing.Days, na.rm=TRUE)) 
# Sum fishing days by year, country, and gear
tab$code <- paste(tab$Country,tab$gear,sep="_") # Create a combined code for country + gear

#### ---- Effort spatial pattern ----
d <- d[d$Geartype %in% unique(gear_table$gear), ] # Filter GFW data to include only gears in gear_table

fp <- d %>% dplyr::group_by(csquares, year, Lat, Lon, Flag, Geartype) %>% dplyr::summarise(effort = sum(Apparent.Fishing.Hours, na.rm=TRUE))
# Calculate total fishing hours by grouping variables

fp <- fp %>% dplyr::group_by(csquares, Lat, Lon, Flag, Geartype) %>% dplyr::summarise(effort = mean(effort, na.rm=TRUE))
# Take the mean effort over the years (if multiple entries exist for each combination)

fp <- fp %>%
  left_join(tab_country, by = c("Flag" = "country_code")) # Merge to get the country name from the code
fp$code <- paste(fp$country,fp$Geartype,sep="_") # Create a combined code
fp <- data.frame(fp) # Convert to data frame

fp_st <- fp %>% 
  dplyr::group_by(code) %>%  # Group by the combined code
  dplyr::mutate(prop = effort / sum(effort, na.rm=TRUE)) %>%  # Calculate the proportion of effort for each code
  dplyr::ungroup() # Ungroup

fp <- as.data.frame(fp_st) # Finalize the data frame

#-----------
fp_codes <- unique(fp$code) # Unique codes in GFW data
fdi_codes <- unique(tab$code) # Unique codes in FDI data
fp_codes # Print to see the codes in fp
fdi_codes # Print to see the codes in fdi
fdi_codes[which(!fdi_codes %in% fp_codes )] # Check which FDI codes are not in fp
fp_codes[which(!fp_codes %in% fdi_codes )] # Check which fp codes are not in FDI

common_codes <- intersect(fp_codes, fdi_codes) # Intersection of codes
print(common_codes) # Print common codes

tab <- tab[tab$code %in% common_codes, ] # Filter FDI table to only the common codes
fp <- fp[fp$code %in% common_codes, ]    # Filter GFW data to only the common codes

# #-----------
# # The following lines (commented out) would create a map with facets by Flag and Geartype
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

y=1 # Initialize a counter (not strictly necessary, but used in the loop)
exp <- list() # Create an empty list to store results

for (y in 1:length(ys)) {
  
  temp <- fp
  colnames(temp)[which(colnames(temp)=="effort")] <- "effort_gfw" # Rename 'effort' to 'effort_gfw'
  taby <- tab[tab$Year == ys[y], c("Year","effort","code")]       # Filter FDI table for the given year
  temp <- temp %>% left_join(taby, by = c("code" = "code"))       # Merge with FDI data
  temp <- temp[!is.na(temp$effort), ]                             # Exclude rows with no FDI effort
  temp$FD <- temp$prop * temp$effort                              # Redistribute the total FDI effort using GFW proportions
  exp[[length(exp)+1]] <- temp                                    # Store the result in the list
  rm(temp)                                                        # Remove the temporary variable
  
}

dd <- do.call(rbind,exp) # Combine all yearly data frames
dd <- data.frame(dd)      # Convert to data frame
head(dd)                  # Inspect the first rows

dd_aggregated <- dd %>% 
  dplyr::group_by(csquares, Lat, Lon, Year) %>% 
  dplyr::summarise(FD_total = sum(FD, na.rm=TRUE)) %>% 
  dplyr::ungroup()
# Summarize total fishing days (FD) per csquare, lat, lon, and year

print(head(dd_aggregated)) # Print a sample of the aggregated data

dd_aggregated$depth <- raster::extract(depth, dd_aggregated[,c("Lon","Lat")]) # Extract the depth from the raster at each point
dd_aggregated$depth  <- -1 * dd_aggregated$depth  # Convert depth values to positive (if negative)
dd_aggregated <- dd_aggregated[dd_aggregated$depth <= 1000, ] # Keep only points with depth <= 1000m
dd_aggregated <- dd_aggregated[!is.na(dd_aggregated$Year) & dd_aggregated$FD_total > 0,]
# Filter out missing years and FD=0

# Aggregate by year to get total FD
dd_aggregated2 <- dd_aggregated %>% 
  group_by(Year) %>% 
  summarise(FD = sum(FD_total, na.rm=TRUE)) %>%
  filter(!is.na(Year))

# Plot a regression of total FD over the years
plot_regression <- ggplot(dd_aggregated2, aes(x = Year, y = FD)) +
  geom_point() +  
  geom_smooth(method = "lm", formula = y ~ x, color = "blue", se = TRUE) + 
  theme_minimal() +
  labs(title = "Linear Regression of FD over Years",
       x = "Year",
       y = "FD (Total)")

print(plot_regression) # Display the plot

cor.test(dd_aggregated2$Year, dd_aggregated2$FD, method = "spearman") # Spearman correlation test

d2013 <- dd_aggregated[dd_aggregated$Year == 2013, ] # Filter rows for 2013
d2012 <- d2013                                       # Create a copy labeled as 2012
d2012$Year <- 2012

dd_aggregated <- rbind(dd_aggregated, d2012) # Add 2012 data (same as 2013) to the dataset
dd_aggregated$fleet <- fleet                 # Tag with fleet type

# Load world map for plotting
world <- ne_countries(scale = "medium", returnclass = "sf") 
# Define Mediterranean extent for mapping
med_extent <- c(xmin = -6, xmax = 37, ymin = 30, ymax = 48) 
label <- paste("Fishing Effort of",fleet) # Title label

# Create a facet_wrap map by year
plot_effort <- ggplot() +
  theme_minimal()+
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

print(plot_effort) # Display the map

# (Plot repeated three times, possibly for iterative checks)
plot_effort <- ggplot() +
  theme_minimal()+
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

print(plot_effort)

# Save the final plot
output_file <- file.path(resdir, paste0("GFW_FDI_rescaled_effort_", paste0(range(ys), collapse="_"), ".jpg"))
ggsave(filename = output_file, plot = plot_effort, width = 10, height = 6, dpi = 300)

#---------------

write.table(dd_aggregated, file.path(resdir, paste0("GFW_FDI_rescaled_effort_",paste0(range(ys),collapse="_"),".csv")),sep=";",row.names=FALSE)
# Export the aggregated data as a CSV

dir.create(file.path(resdir,"rasters")) # Create a directory for storing rasters

# Loop over unique years in dd_aggregated to create yearly raster
for (year in unique(dd_aggregated$Year)) {
  
  dd_year <- dd_aggregated %>% filter(Year == year) # Filter data for that year
  
  # Create a raster from Lon, Lat, and FD_total
  r <- rasterFromXYZ(dd_year[, c("Lon", "Lat", "FD_total")])
  
  # Set the coordinate reference system to WGS84
  crs(r) <- CRS("+proj=longlat +datum=WGS84")
  
  # Build output raster filename
  output_file <- file.path(resdir,"rasters", paste0(year,"_FE_",fleet, ".tif"))
  
  # Save raster
  writeRaster(r, output_file, format = "GTiff", overwrite = TRUE)
  
  message(paste("Raster saved for year:", year)) # Console message
}
