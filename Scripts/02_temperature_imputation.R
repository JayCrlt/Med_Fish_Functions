#### Setting up          ----
library(dplyr) ; library(sf) ; library(lubridate) ; library(FNN) ; library(pbapply) ; library(parallel)

## Functions
source("Scripts/00_functions_script.R")

## Download data
sp_code_list       <- read.delim("Data/MEDITS_spp.codes.csv", sep = ";")
TATB_WMED          <- read.delim("Data/TATB_WMED_1999-2021_clean.csv", sep = ";")
TATB_EMED          <- read.delim("Data/TATB_EMED_1999-2021_clean.csv", sep = ";")
Medits_total       <- rbind(TATB_WMED, TATB_EMED) |> 
  left_join(Temp_MED <- read.delim("Data/HAULS_DB_ENV_HINDCAST_COMPLETE_2024.csv", sep = ";") |> 
              select(-id) |> rename(id = ids) |> select(c(id,botTemp_anom)) |> rename(Temp_anomaly = botTemp_anom))
Hexagonal_grid     <- st_read("Data/Grid_0-1000m_Med.shp")

# add spatial hex info
medits_sf               <- st_as_sf(data.frame(Longitude = as.numeric(gsub(",", ".", Medits_total$MEAN_LONGITUDE_DEC)), 
                                               Latitude = as.numeric(gsub(",", ".", Medits_total$MEAN_LATITUDE_DEC))), 
                                    coords = c("Longitude", "Latitude"), crs = 4326)
medits_sf_hex           <- st_join(medits_sf, Hexagonal_grid["grid_id"], left = TRUE, join = st_nearest_feature)
Medits_total$HEX_ID     <- medits_sf_hex$grid_id
Medits_total$MASS_IND   <- Medits_total$TOTAL_WEIGHT_IN_THE_HAUL / Medits_total$TOTAL_NUMBER_IN_THE_HAUL

# Change Temperature Variable
MEDIT_TEMP       <- cbind(Medits_total, medits_sf) |> data.frame() |> 
  mutate(TEMP = rowMeans(across(c(BOTTOM_TEMPERATURE_BEGINNING, BOTTOM_TEMPERATURE_END), ~ {
    temp_num     <- as.numeric(gsub(",", ".", .)); ifelse(temp_num > 100, temp_num / 10, temp_num)}), na.rm = TRUE))
MEDIT_TEMP$Date  <- as.Date(paste(MEDIT_TEMP$YEAR, MEDIT_TEMP$MONTH, "01", sep = "-"))
MEDIT_TEMP$Month <- MEDIT_TEMP$MONTH
temp_available   <- MEDIT_TEMP |> filter(!is.na(TEMP)) |> select(HEX_ID, Date, Month, TEMP)
hex_centroids    <- st_centroid(Hexagonal_grid)
hex_coords       <- st_coordinates(hex_centroids)
hex_lookup       <- data.frame(HEX_ID = Hexagonal_grid$grid_id, X = hex_coords[,1], Y = hex_coords[,2])

# Apply to dataset
MEDIT_TEMP$TEMP_imputed <- pbapply::pbmapply(impute_temp, MEDIT_TEMP$HEX_ID, MEDIT_TEMP$Date, MEDIT_TEMP$Month)
MEDIT_TEMP = MEDIT_TEMP |> mutate(TEMP = TEMP_imputed) |> select(-c(geometry, Date, Month, TEMP_imputed))

# Save
save(MEDIT_TEMP, file = "Outputs/dat_proc/Medit_Temp.RData")