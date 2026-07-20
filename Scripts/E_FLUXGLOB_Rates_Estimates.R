#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("rphylopic") ; library("FNN") ; library("pbapply")

## Charge from previous scripts
load("Outputs/FLUXGLOB/dat_proc/FLUXGLOB_with_temp.RData")
load("Outputs/FLUXGLOB/dat_proc/All_Species_FLUXGLOB.RData")
## Functions
source("Scripts/00_functions_script.R")

# Work clean, remove haul_ID with no num and no weight
FLUXGLOB_cor = FLUXGLOB_cor |> mutate(Haul_swept_area = clean_numeric(Haul_swept_area))
FLUXGLOB_cor_clean <- FLUXGLOB_cor |> filter(!is.na(Num_CPUA), !is.na(Weight_CPUA), !(Num_CPUA == 0 & Weight_CPUA == 0)) |>
  group_by(Haul_ID) |> filter(!any(xor(Num_CPUA == 0, Weight_CPUA == 0))) |> ungroup()

# Impute Haul swept area and duration with mean observation
hauls <- FLUXGLOB_cor_clean |> distinct(Haul_ID, .keep_all = T) |>
  mutate(Haul_swept_area = clean_numeric(Haul_swept_area), Haul_duration   = clean_numeric(Haul_duration)) |>
  group_by(Dataset, PROVINCE, Year) |>
  mutate(Haul_swept_area = if_else(is.na(Haul_swept_area), mean(Haul_swept_area, na.rm = T), Haul_swept_area),
    Haul_duration = if_else(is.na(Haul_duration), mean(Haul_duration, na.rm = T), Haul_duration)) |> ungroup() |>
  group_by(Dataset, PROVINCE) |>
  mutate(Haul_swept_area = if_else(is.na(Haul_swept_area), mean(Haul_swept_area, na.rm = T), Haul_swept_area),
    Haul_duration = if_else(is.na(Haul_duration), mean(Haul_duration, na.rm = T), Haul_duration)) |> ungroup() |>
  group_by(Dataset) |>
  mutate(Haul_swept_area = if_else(is.na(Haul_swept_area), mean(Haul_swept_area, na.rm = T), Haul_swept_area),
    Haul_duration = if_else(is.na(Haul_duration), mean(Haul_duration, na.rm = T), Haul_duration)) |> ungroup()
FLUXGLOB_cor_clean <- FLUXGLOB_cor_clean |> dplyr::select(-Haul_swept_area, -Haul_duration) |>
  left_join(st_drop_geometry(hauls |> dplyr::select(Haul_ID, Haul_swept_area, Haul_duration)), by = "Haul_ID")

# Some species are unexpectedly heavy...
FLUXGLOB_cor_clean = FLUXGLOB_cor_clean |> mutate(Num = Num_CPUA * Haul_swept_area, 
                                      Weight = Weight_CPUA * Haul_swept_area,
                                      Fish_weight_avg = Weight / Num)

### Let's check total biomass
Species                  <- gsub(" ", "_", tab_Glob$Species)
phy_weight               <- fishtree_phylogeny(species = unique(Species))
phy_weight$tip.label     <- gsub(" ", "_", phy_weight$tip.label)

tab_Glob_weight = tab_Glob |> 
  left_join(rfishbase::load_taxa() |> dplyr::select(SpecCode, Species), by = "Species") |> 
  mutate(SpecCode = as.character(SpecCode)) |> 
  left_join((rfishbase::popchar() |> dplyr::select(SpecCode, Wmax)) |> mutate(SpecCode = as.character(SpecCode)), 
            by = "SpecCode", relationship = "many-to-many") |> 
  group_by(Species) |> 
  slice_max(order_by = coalesce(Wmax, -Inf), n = 1, with_ties = FALSE) |> ungroup() |> 
  left_join((rfishbase::popgrowth() |> dplyr::select(SpecCode, Winfinity)) |> mutate(SpecCode = as.character(SpecCode)), 
            by = "SpecCode", relationship = "many-to-many") |> 
  group_by(Species) |> slice_max(order_by = coalesce(Winfinity, -Inf), n = 1, with_ties = FALSE) |> ungroup() |> 
  impute_trait_phylopars_FG(trait_col = "Wmax", phy = phy_weight) |> fill_trait_hierarchy_FG("Wmax") |> 
  impute_trait_phylopars_FG(trait_col = "Winfinity", phy = phy_weight) |> fill_trait_hierarchy_FG("Winfinity") |> 
  dplyr::select(c(Species, Wmax, Winfinity)) |> 
  right_join(tab_Glob, by = "Species")

### Merge both datasets
trait_cols <- tab_Glob_weight |> dplyr::select(Species_FLUXGLOB, Qc, Qn, Qp, Dc, Dn, Dp, h, ac, an, ap,
    F0Nz, F0Pz, mdw, lwa, lwb, linf, K, t0, r, alpha, Wmax, Winfinity)
FLUXGLOB_cor_clean <- FLUXGLOB_cor_clean |> left_join(trait_cols, by = "Species_FLUXGLOB") |> 
  dplyr::select(-c(subplot_Fig1, Class, SuperClass)) |>
  mutate(SBTemp = clean_numeric(SBTemp), SBTemp = if_else(SBTemp > 40, NA_real_, SBTemp), 
         SBTemp = if_else(SBTemp <= 0, NA_real_, SBTemp))

### Impute missing observed temperature
hauls                 <- FLUXGLOB_cor_clean |> distinct(Haul_ID, .keep_all = TRUE)
hauls$Month           <- as.integer(hauls$Month)
coords                <- sf::st_coordinates(hauls)
temp_available        <- hauls |> dplyr::filter(!is.na(SBTemp))
temp_xy               <- sf::st_coordinates(temp_available)
temp_available$X      <- temp_xy[, 1]
temp_available$Y      <- temp_xy[, 2]
missing               <- which(is.na(hauls$SBTemp))
hauls$SBTemp[missing] <- pbsapply(missing, function(i) impute_temp_FG(hauls$Month[i], hauls$Year[i], coords[i, ]))
FLUXGLOB_cor_clean    <- FLUXGLOB_cor_clean |> 
  left_join(hauls |> st_drop_geometry() |> dplyr::select(Haul_ID, SBTemp), by = "Haul_ID", suffix = c("", ".new")) |> 
  mutate(SBTemp = coalesce(SBTemp.new, SBTemp)) |> dplyr::select(-SBTemp.new)

# Clean errors according to maximum weight
FLUXGLOB_final <- FLUXGLOB_cor_clean |> 
  mutate(NegativeDifference = (if_else(Winfinity > Wmax, Winfinity, Wmax) - Fish_weight_avg) < -100) |>
  group_by(Haul_ID) |> dplyr::filter(!any(NegativeDifference, na.rm = TRUE)) |> ungroup() |> # Remove haul ID when 1 obs failed
  #dplyr::filter(NegativeDifference == "FALSE") |>                                           # Remove only obs failing
  mutate(SMR = 0.002 * Fish_weight_avg^0.836,
         MMR = 0.006 * Fish_weight_avg^0.779,
         f0 = SMR / Fish_weight_avg^0.836 *
           exp(-0.138 * (1 / (botTemp + 273.15) - 1 / (botTemp + (botTemp - SBTemp) + 273.15)) / 8.617e-5),
         f0 = if_else(is.nan(f0), mean(f0[!is.nan(f0)], na.rm = TRUE), f0),
         theta = (SMR + MMR) / (2 * SMR),
         Haul_swept_area = clean_numeric(Haul_swept_area)) |> relocate(mdw, .after = lwb) |>
  mutate(TL = (Fish_weight_avg / lwa)^(1 / lwb)) |> dplyr::select(-NegativeDifference)

### Not enough taxonomic details – Remove
FLUXGLOB_final <- FLUXGLOB_final |> filter(!Haul_ID %in%
           (FLUXGLOB_final |> filter(is.na(Wmax) | is.na(Winfinity)) |> pull(Haul_ID))) |>
  mutate(across(where(is.numeric), unname))

## Define data reduction  
sprintf("Data reduction of %.2f%%", 100 * (1 - nrow(FLUXGLOB_final) / nrow(FLUXGLOB_cor)))

## FLUXGLOB loop workflow
FLUXGLOB_final = FLUXGLOB_final |> 
  mutate(Fn_mean = NA_real_, Fn_Q1   = NA_real_, Fn_Q3   = NA_real_, Fp_mean = NA_real_, 
         Fp_Q1   = NA_real_, Fp_Q3   = NA_real_, Gc_mean = NA_real_, Gc_Q1   = NA_real_, 
         Gc_Q3   = NA_real_, Ic_mean = NA_real_, Ic_Q1   = NA_real_, Ic_Q3   = NA_real_)
ids <- sort(unique(FLUXGLOB_final$ID_Fig1))
for(id in ids){
  cat("\nProcessing ID", id, "\n")
  idx <- which(FLUXGLOB_final$ID_Fig1 == id)
  dat <- FLUXGLOB_final[idx, ]
  pb  <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
  for(i in seq_len(nrow(dat))){
    param <- list(
      ac_m    = dat$ac[i],
      ac_sd   = dat$ac[i] * 0.05,
      an_m    = dat$an[i],
      an_sd   = dat$an[i] * 0.05,
      ap_m    = dat$ap[i],
      ap_sd   = dat$ap[i] * 0.05,
      r_m     = dat$r[i],
      f0_m    = dat$f0[i],
      theta_m = dat$theta[i],
      Dc_m    = dat$Dc[i],
      Dc_sd   = dat$Dc[i] * 0.05,
      Dn_m    = dat$Dn[i],
      Dn_sd   = dat$Dn[i] * 0.05,
      Dp_m    = dat$Dp[i],
      Dp_sd   = dat$Dp[i] * 0.05,
      k_m     = dat$K[i],
      linf_m  = dat$linf[i],
      lwa_m   = dat$lwa[i],
      lwb_m   = dat$lwb[i],
      Qc_m    = dat$Qc[i],
      Qc_sd   = dat$Qc[i] * 0.05,
      Qn_m    = dat$Qn[i],
      Qn_sd   = dat$Qn[i] * 0.05,
      Qp_m    = dat$Qp[i],
      Qp_sd   = dat$Qp[i] * 0.05,
      t0_m    = dat$t0[i],
      h_m     = dat$h[i],
      F0nz_m  = dat$F0Nz[i],
      F0pz_m  = dat$F0Pz[i],
      mdw_m   = dat$mdw[i],
      v_m     = dat$SBTemp[i],
      alpha_m = 0.836)
    tryCatch({invisible(capture.output(
          res <- extract(cnp_model_mcmc(TL = dat$TL[i], param = param), c("Fn", "Fp", "Gc", "Ic"))))
      dat$Fn_mean[i] <- res$Fn_mean
      dat$Fn_Q1[i]   <- res$`Fn_2.5%`
      dat$Fn_Q3[i]   <- res$`Fn_97.5%`
      dat$Fp_mean[i] <- res$Fp_mean
      dat$Fp_Q1[i]   <- res$`Fp_2.5%`
      dat$Fp_Q3[i]   <- res$`Fp_97.5%`
      dat$Gc_mean[i] <- res$Gc_mean
      dat$Gc_Q1[i]   <- res$`Gc_2.5%`
      dat$Gc_Q3[i]   <- res$`Gc_97.5%`
      dat$Ic_mean[i] <- res$Ic_mean
      dat$Ic_Q1[i]   <- res$`Ic_2.5%`
      dat$Ic_Q3[i]   <- res$`Ic_97.5%`}, error = function(e){
      message("Error in row ", i, " (ID ", id, "): ", conditionMessage(e))})
    setTxtProgressBar(pb, i)}
  close(pb)
  cols <- c("Fn_mean","Fn_Q1","Fn_Q3", "Fp_mean","Fp_Q1","Fp_Q3", "Gc_mean","Gc_Q1","Gc_Q3", "Ic_mean","Ic_Q1","Ic_Q3")
  for(cl in cols){FLUXGLOB_final[[cl]][idx] <- dat[[cl]]}
  # Save checkpoint after finishing this ID
  save(FLUXGLOB_final, file = paste0("Outputs/FLUXGLOB/dat_proc/FLUXGLOB_checkpoint_ID", id, ".RData"))
  rm(dat)
  gc()}

#### Export the data  ----
save(FLUXGLOB_final, file = "Outputs/FLUXGLOB/dat_proc/FLUXGLOB_with_fluxes.RData")