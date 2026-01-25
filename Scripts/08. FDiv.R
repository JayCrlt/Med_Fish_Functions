#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom") ; library("mFD") ; library("funspace")

## Download data
sp_Traits <- readxl::read_xlsx("Data/Traits_MEDITS_Fun.xlsx") 
traits_df <- sp_Traits |> dplyr::select(SPECIES, lengthC, lifespanC, vertical, diet, TempPref)
traits_mat <- traits_df |> tibble::column_to_rownames("SPECIES") |> 
  mutate(lengthC = ordered(lengthC), lifespanC = ordered(lifespanC), TempPref = ordered(TempPref), 
         vertical = factor(vertical), diet = factor(diet))
trait_types <- data.frame(trait_name = colnames(traits_mat), trait_type = c("O", "O", "N", "N", "O"))
dist_funct  <- mFD::funct.dist(sp_tr = traits_mat, tr_cat = trait_types, metric = "gower", scale_euclid = "none")
funct_space <- mFD::quality.fspaces(sp_dist = dist_funct, maxdim_pcoa = 3, deviation_weighting = "absolute", fdist_scaling = FALSE)
funct_space$details_fspaces

sp_coords <- funct_space$details_fspaces$sp_pc_coord
fs_coords <- sp_coords[, 1:2]
fs_coords_mat <- as.matrix(fs_coords)
fs_object <- funspace::funspace(fs_coords_mat)
fs_object <- funspace::funspace(fs_coords_mat)
plot(fs_object, type = "global", pnt = T, pnt.col = "black", pnt.cex = 0.7, quant.plot = T, quant = c(0.5, 0.95, 0.999), 
     axis.title.x = "PCoA 1", axis.title.y = "PCoA 2")

fs_coords_mat_2_3 <- as.matrix(fs_coords[, 2:3])  
fs_object_2_3 <- funspace::funspace(fs_coords_mat_2_3)
plot(fs_object_2_3, type = "global", pnt = T, pnt.col = "black", pnt.cex = 0.7, quant.plot = TRUE, quant = c(0.5, 0.95, 0.999),
  axis.title.x = "PCoA 2", axis.title.y = "PCoA 3")