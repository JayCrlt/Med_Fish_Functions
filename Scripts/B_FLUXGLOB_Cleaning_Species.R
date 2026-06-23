#### Setting up          ----
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales") 
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("ggforce")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("purrr")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("hexbin") ; library("leaflet.extras") 
library("RColorBrewer") ; library("MASS") ; library("sf") ; library("brms") ; library("stringr") ; library("ggtext")

## Functions
source("Scripts/00_functions_script.R")
source("Scripts/A_FLUXGLOB_Context.R")

## Clean Species name
Species_list       <- FLUXGLOB |> distinct(Species) |> left_join(rfishbase::load_taxa())
Species_list_NA_FB <- Species_list |> dplyr::filter(is.na(SpecCode)) |> pull(Species)
### 1) Wrong denomination due to out-to-date names
Species_list_cor   <- data.frame(Species_FLUXGLOB = Species_list_NA_FB, Species = validate_names(Species_list_NA_FB)) |>
  left_join(rfishbase::load_taxa()) |> mutate(Species_FLUXGLOB = case_when(Species_FLUXGLOB == "Trigloporus" ~ "Chelidonichthys",
                                                                           Species_FLUXGLOB == "Brotulidae" ~ "Brotulinae",
                                                                           Species_FLUXGLOB == "Percoidei" ~ "Perciformes/Percoidei",
                                                                           TRUE ~ Species_FLUXGLOB),
                                              Class = case_when(Species_FLUXGLOB == "Selachii" ~ "Elasmobranchii", TRUE ~ Class),
                                              n_words = stringr::str_count(Species_FLUXGLOB, "\\S+")) |> 
  distinct(Species_FLUXGLOB, Species, SpecCode, Genus, Subfamily, Family, Order, Class, SuperClass, n_words)
### 2) ID to genus, Subfamily, Family, Order level resolution
taxa               <- load_taxa()
order_lookup       <- taxa |> distinct(Order, Class)
family_lookup      <- taxa |> distinct(Family, Order, Class)
subfamily_lookup   <- taxa |> distinct(Subfamily, Family, Order, Class)
genus_lookup       <- taxa |> distinct(Genus, Subfamily, Family, Order, Class)
remaining_order    <- Species_list_cor |> dplyr::filter(n_words == 1) |> dplyr::select(-c(Class)) |> 
  left_join(order_lookup, by = c("Species_FLUXGLOB" = "Order")) |> mutate(Species = NA_character_) |> 
  dplyr::select(names(Species_list_cor)) |> drop_na(Class) |> mutate(Order = Species_FLUXGLOB)
remaining_family   <- Species_list_cor |> dplyr::filter(n_words == 1) |> dplyr::select(-c(Order, Class)) |> 
  left_join(family_lookup, by = c("Species_FLUXGLOB" = "Family")) |> mutate(Species = NA_character_) |> 
  dplyr::select(names(Species_list_cor)) |> drop_na(Class) |> mutate(Family = Species_FLUXGLOB)
remaining_subfam   <- Species_list_cor |> dplyr::filter(n_words == 1) |> dplyr::select(-c(Order, Class, Family)) |> 
  left_join(subfamily_lookup, by = c("Species_FLUXGLOB" = "Subfamily")) |> mutate(Species = NA_character_) |> 
  dplyr::select(names(Species_list_cor)) |> drop_na(Class) |> mutate(Subfamily = Species_FLUXGLOB)
remaining_genus    <- Species_list_cor |> dplyr::filter(n_words == 1) |> dplyr::select(-c(Order, Class, Family, Subfamily)) |> 
  left_join(genus_lookup, by = c("Species_FLUXGLOB" = "Genus")) |> mutate(Species = NA_character_) |> 
  dplyr::select(names(Species_list_cor)) |> drop_na(Class) |> mutate(Genus = Species_FLUXGLOB)
remaining = rbind(remaining_order, remaining_family, remaining_subfam, remaining_genus) 
Species_list_nom   <- Species_list_cor |> drop_na(Class)
Species_list_nom   <- rbind(Species_list_nom, remaining)
### 3) Genus spp.
Species_list_gen   <- Species_list_cor |> dplyr::filter(Species_FLUXGLOB %in% 
  Species_list_NA_FB[grepl("\\s+spp?\\.?$", Species_list_NA_FB)]) |>
  dplyr::select(Species_FLUXGLOB) |> mutate(Genus = sub("\\s+spp?\\.?$", "", Species_FLUXGLOB)) |> 
  left_join(genus_lookup, by = "Genus") |> mutate(Species = NA_character_, SpecCode = NA_character_, SuperClass = NA_character_,
  n_words = 2) |> dplyr::select(names(Species_list_cor))
Species_list_cor   <- rbind(Species_list_nom, Species_list_gen) |> dplyr::select(-n_words)
### Non-identified species / Manual data wrangling
Species_list_OSNI  <- data.frame(Species_FLUXGLOB = c("Symphodus doderleini ", "Allinectes ectenes", "Hyporthodus drummondhayi", 
                                                      "Hypocritichthys analis", "Sinosturio medirostris", "Bothrocara zestum"), 
                                 Species = c("Symphodus doderleini", "Careproctus ectenes", "Epinephelus drummondhayi", 
                                             "Hyperprosopon anale", "Acipenser medirostris", "Bothrocara brunneum")) |> 
  left_join(rfishbase::load_taxa()) |> dplyr::select(names(Species_list_cor))
Species_list_cor   <- rbind(Species_list_cor, Species_list_OSNI) 
Species_list_cor$Species_FLUXGLOB[Species_list_cor$Species_FLUXGLOB == "Perciformes/Percoidei"] = "Percoidei"
Species_list_cor$Species_FLUXGLOB[Species_list_cor$Species_FLUXGLOB == "Brotulinae"] = "Brotulidae"
dupli              <- Species_list_cor[Species_list_cor$Species_FLUXGLOB == "Chelidonichthys",] 
dupli$Species_FLUXGLOB[dupli$Species_FLUXGLOB == "Chelidonichthys"] = "Trigloporus"
Species_list_cor   <- rbind(Species_list_cor, dupli) 
(OSNI = setdiff(Species_list_NA_FB, Species_list_cor$Species_FLUXGLOB)) ; length(OSNI)
### Only 5 individuals named as "Perciformes" in the total dataset – Will have to be removed
Species_list_no_NA <- Species_list |> mutate(Species_FLUXGLOB = Species) |> dplyr::select(names(Species_list_cor)) |> 
  drop_na(SpecCode)
Species_list       <- rbind(Species_list_no_NA, Species_list_cor) |> arrange(Species_FLUXGLOB) 

### Final Clean FLUXGLOB dataset
FLUXGLOB_cor = FLUXGLOB |> dplyr::select(-c(SpecCode, Family, Genus)) |> 
  left_join(Species_list |> dplyr::filter(!(Species_FLUXGLOB == "Leptocephalus" & is.na(Subfamily) & Family == "Notacanthidae")), 
            by = c("Species" = "Species_FLUXGLOB")) |> dplyr::select(-Species) |> rename(Species = Species.y)

##### Save NMDS Model
save(FLUXGLOB_cor, file = "Outputs/FLUXGLOB/dat_proc/FLUXGLOB.RData")
save(Species_list, file = "Outputs/FLUXGLOB/dat_proc/FLUXGLOB_Species_list.RData")