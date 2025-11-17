#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales") ; library("sf") 
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("bayestestR")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom") ; library("gstat") ; library("forcats")

## Download data
Env <- read.delim("Data/HAULS_DB_ENV_HINDCAST_COMPLETE_2024.csv", sep = ";")
FD  <- read.delim("Data/Grid_taxfd.csv", sep = ",")
reg <- st_read("Data/reg_/reg_.shp") |> mutate(Name = coalesce(Name, layer)) |> select(Name, geometry)

## Charge from previous scripts
load("Data/HAULS_DB_FPI.RData")
load("Outputs/dat_proc/medits_sf_percentile.Rdata")
load("Outputs/dat_proc/Medit_FunCatch_without_NA.RData")
load("Outputs/dat_proc/Temp_change_slopes.RData")

## Functions
source("Scripts/00_functions_script.R")

Env_FI = merge(Env, HAULS_DB_FPI, by = c("ids", "X", "Y", "YEAR")) |> 
  mutate(GSA = as.character(GSA), MONTH = as.character(MONTH))
merged_sf <- st_join(Slopes, st_as_sf(Env_FI, coords = c("X", "Y"), crs = 4326), join = st_nearest_feature) |> 
  dplyr::filter(Mf_trend_class %in% c("4. significantly positive", "1. significantly negative")) |> 
  mutate(FPI_combined = coalesce(nearest_FPI, FPI_tot)) |> select(-c(nearest_FPI, FPI_tot)) |> 
  mutate(Mf_bin = ifelse(Multifunctionality_slope > 0, 1, 0), abs_botTemp_anom = abs(botTemp_anom),
         biom_bin = ifelse(Biomass_slope > 0, 1, 0)) |> rename(grid.id = HEX_ID) |>
  left_join((FD |> pivot_longer(cols = c(fdiv, feve, fric, fvuln, shannon_diversity, sp_richness), names_to = "metric", values_to = "value") |>
               group_by(grid.id, metric) |> filter(!is.na(value) & !is.na(YEAR)) |> filter(n() >= 2) |> do(tidy(lm(value ~ YEAR, data = .))) |>
               filter(term == "YEAR") |> select(grid.id, metric, slope = estimate)) |> 
              pivot_wider(names_from = metric, values_from = slope, names_glue = "{metric}_slope"), by = "grid.id") |> 
  left_join(Medit_FunCatch_without_NA |> pivot_longer(cols = c(TOTAL_NUMBER_IN_THE_HAUL), names_to = "metric", values_to = "value") |>
              group_by(HEX_ID, metric) |> filter(!is.na(value) & !is.na(YEAR)) |> filter(n() >= 2) |> do(tidy(lm(value ~ YEAR, data = .))) |>
              filter(term == "YEAR") |> select(HEX_ID, metric, slope = estimate) |> rename("grid.id" = "HEX_ID") |> 
              pivot_wider(names_from = metric, values_from = slope, names_glue = "{metric}_slope"), by = "grid.id") |> st_as_sf() |> 
  rename("abundance_change" = "TOTAL_NUMBER_IN_THE_HAUL_slope", "fdiv_change" = "fdiv_slope")

## Asure that regions works
# st_is_valid(reg) ; reg <- st_make_valid(reg)
# merged_with_regions <- st_join(merged_sf, reg, join = st_within)

## Check correlation before modelling
merged_sf |> st_drop_geometry() |> select(c(2:9, 65:79)) |> mutate(across(everything(), as.numeric)) |> 
  cor(use = "pairwise.complete.obs") |> corrplot::corrplot(method = "color", type = "upper",  addCoef.col = "black", tl.cex = .7, number.cex = .5)

## General Model
model = brms::brm(formula = Mf_bin ~ botTemp_anom + DEPTH + chl + FPI_combined + fdiv_change + (1|GSA), 
                  data = merged_sf, cores = 3, chains = 3, iter = 10000, family = bernoulli(link = "logit"))
(post_summary <- as.data.frame(fixef(model)))
bayes_R2(model)

## Define the contribution of each variable
model_std <- brm(Mf_bin ~ scale(botTemp_anom) + scale(DEPTH) + scale(chl) + scale(FPI_combined) + scale(fdiv_change) + (1|GSA),
                 data = merged_sf, family = bernoulli(link = "logit"), cores = 3, chains = 3, iter = 10000)
fixef(model_std)

draws_df <- as_draws_df(model_std) |> select(matches("^b_")) |> pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") |>
  mutate(parameter = gsub("^b_", "", parameter)) |> group_by(parameter) |>
  summarise(median = median(value), hdi_50_low  = hdi(value, ci = 0.50)$CI_low, hdi_50_high = hdi(value, ci = 0.50)$CI_high,
    hdi_95_low  = hdi(value, ci = 0.925)$CI_low, hdi_95_high = hdi(value, ci = 0.925)$CI_high, .groups = "drop") |> 
  mutate(parameter = recode(parameter, "b_Intercept" = "Intercept", "scaleDEPTH" = "Depth", "scaleFPI_combined" = "Fishing Pressure Index",
                            "scalebotTemp_anom" = "Sea Bottom Temperature Anomaly", "scalechl" = "Chlorophyll-a concentration",
                            "scalefdiv_change" = "Functional diversity change")) |> 
  mutate(fill_color = c("#ffffff","#FF6666","#ffffff","#ffffff","#9999FF","#9999FF"),
         seg_color = c("#9999FF","#FF6666","#FF6666","#9999FF","#9999FF","#9999FF"))

Figure_S5 <- ggplot(draws_df, aes(x = median, y = fct_reorder(parameter, median))) +
  geom_segment(aes(x = hdi_95_low, xend = hdi_95_high, yend = parameter, color = seg_color), size = 0.5) +
  geom_segment(aes(x = hdi_50_low, xend = hdi_50_high, yend = parameter, color = seg_color), size = 2) +
  geom_point(aes(fill = fill_color), shape = 21, size = 5, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  theme_classic() + scale_fill_identity() + scale_color_identity() + labs(x = "Effect size", y = "") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text       = element_text(size = 16),
       legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")

## Conditional effect Fishing Pressure
ce_FPI <- conditional_effects(model, effects = "FPI_combined")
ce_FPI$FPI_combined$gam_low = predict(mgcv::gam(lower__ ~ s(FPI_combined, k = 10), data = ce_FPI$FPI_combined, family = gaussian()), 
                                      newdata = ce_FPI$FPI_combined)
ce_FPI$FPI_combined$gam_upp = predict(mgcv::gam(upper__ ~ s(FPI_combined, k = 10), data = ce_FPI$FPI_combined, family = gaussian()), 
                                      newdata = ce_FPI$FPI_combined)
ce_FPI$FPI_combined$gam_est = predict(mgcv::gam(estimate__ ~ s(FPI_combined, k = 10), data = ce_FPI$FPI_combined, family = gaussian()), 
                                      newdata = ce_FPI$FPI_combined)
poly_low <- ce_FPI$FPI_combined %>% select(FPI_combined, gam_low) %>% rename(y = gam_low) %>% arrange(FPI_combined) %>%
  mutate(ymin = 0) %>% select(FPI_combined, y, ymin)
poly_upp <- ce_FPI$FPI_combined %>% select(FPI_combined, gam_upp) %>% rename(y = gam_upp) %>% arrange(FPI_combined) %>%
  mutate(ymax = 1) %>% select(FPI_combined, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_FPI$FPI_combined$FPI_combined)+0.001, max(ce_FPI$FPI_combined$FPI_combined)-0.001, length.out = 500), 
                             y = seq(0.01, 0.99, length.out = 500)) |> mutate(fill_value = y)

Figure_5A <- ggplot(ce_FPI$FPI_combined, aes(x = FPI_combined)) +
    geom_raster(data = gradient_data, aes(x = x, y = y, fill = fill_value), interpolate = TRUE) +
    scale_fill_gradient(low = "#FF6666", high = "#9999FF") +
    geom_ribbon(data = poly_low, aes(ymin = ymin, ymax = y), fill = "white") +
    geom_ribbon(data = poly_upp, aes(ymin = y, ymax = ymax), fill = "white") +
    geom_line(aes(y = gam_est), size = 1, color = "black") +
    geom_line(aes(y = gam_low), size = 0.5, color = "black") +
    geom_line(aes(y = gam_upp), size = 0.5, color = "black") +
    geom_point(data = merged_sf, aes(y = Mf_bin, x = FPI_combined, fill = Mf_bin), shape = 21, size = 3) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(-1, 0, 1)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
    annotate("text", x = min(ce_FPI$FPI_combined$FPI_combined)+0.01, y = 0.52, 
             label = "no change", hjust = 0, vjust = 0, size = 4, fontface = "italic") +
    coord_cartesian(ylim = c(0, 1)) + theme_classic() +
    labs(x = "Fish Pressure Index", y = "Multifunctional Index change") +
     theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
           plot.title      = element_text(size = 20),
           axis.title      = element_text(size = 18),
           axis.text       = element_text(size = 16),
           legend.title    = element_text(size = 18),
           legend.text     = element_text(size = 16),
           legend.position = "none")  

## Conditional effect Climate Change
ce_sbta <- conditional_effects(model, effects = "botTemp_anom") 
ce_sbta$botTemp_anom$gam_low = predict(mgcv::gam(lower__ ~ s(botTemp_anom, k = 10), data = ce_sbta$botTemp_anom, family = gaussian()), 
                                      newdata = ce_sbta$botTemp_anom)
ce_sbta$botTemp_anom$gam_upp = predict(mgcv::gam(upper__ ~ s(botTemp_anom, k = 10), data = ce_sbta$botTemp_anom, family = gaussian()), 
                                      newdata = ce_sbta$botTemp_anom)
ce_sbta$botTemp_anom$gam_est = predict(mgcv::gam(estimate__ ~ s(botTemp_anom, k = 10), data = ce_sbta$botTemp_anom, family = gaussian()), 
                                      newdata = ce_sbta$botTemp_anom)
poly_low <- ce_sbta$botTemp_anom %>% select(botTemp_anom, gam_low) %>% rename(y = gam_low) %>% arrange(botTemp_anom) %>%
  mutate(ymin = 0) %>% select(botTemp_anom, y, ymin)
poly_upp <- ce_sbta$botTemp_anom %>% select(botTemp_anom, gam_upp) %>% rename(y = gam_upp) %>% arrange(botTemp_anom) %>%
  mutate(ymax = 1) %>% select(botTemp_anom, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_sbta$botTemp_anom$botTemp_anom)+0.001, max(ce_sbta$botTemp_anom$botTemp_anom)-0.001, length.out = 500), 
                             y = seq(0.01, 0.995, length.out = 500)) |> mutate(fill_value = y)

Figure_5B <- ggplot(ce_sbta$botTemp_anom, aes(x = botTemp_anom)) +
  geom_raster(data = gradient_data, aes(x = x, y = y, fill = fill_value), interpolate = TRUE) +
  scale_fill_gradient(low = "#FF6666", high = "#9999FF") +
  geom_ribbon(data = poly_low, aes(ymin = ymin, ymax = y), fill = "white") +
  geom_ribbon(data = poly_upp, aes(ymin = y, ymax = ymax), fill = "white") +
  geom_line(aes(y = gam_est), size = 1, color = "black") +
  geom_line(aes(y = gam_low), size = 0.5, color = "black") +
  geom_line(aes(y = gam_upp), size = 0.5, color = "black") +
  geom_point(data = merged_sf, aes(y = Mf_bin, x = botTemp_anom, fill = Mf_bin), shape = 21, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + theme_classic() +
  labs(x = "Bottom Sea Temperature Anomaly", y = "") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text.x     = element_text(size = 16),
        axis.text.y     = element_blank(),      
        axis.ticks.y    = element_blank(),    
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")  

## Conditional effect depth
ce_depth <- conditional_effects(model, effects = "DEPTH") 
ce_depth$DEPTH$gam_low = predict(mgcv::gam(lower__ ~ s(DEPTH, k = 10), data = ce_depth$DEPTH, family = gaussian()), 
                                       newdata = ce_depth$DEPTH)
ce_depth$DEPTH$gam_upp = predict(mgcv::gam(upper__ ~ s(DEPTH, k = 10), data = ce_depth$DEPTH, family = gaussian()), 
                                       newdata = ce_depth$DEPTH)
ce_depth$DEPTH$gam_est = predict(mgcv::gam(estimate__ ~ s(DEPTH, k = 10), data = ce_depth$DEPTH, family = gaussian()), 
                                       newdata = ce_depth$DEPTH)
poly_low <- ce_depth$DEPTH %>% select(DEPTH, gam_low) %>% rename(y = gam_low) %>% arrange(DEPTH) %>%
  mutate(ymin = 0) %>% select(DEPTH, y, ymin)
poly_upp <- ce_depth$DEPTH %>% select(DEPTH, gam_upp) %>% rename(y = gam_upp) %>% arrange(DEPTH) %>%
  mutate(ymax = 1) %>% select(DEPTH, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_depth$DEPTH$DEPTH)+0.001, max(ce_depth$DEPTH$DEPTH)-0.001, length.out = 500), 
                             y = seq(0.01, 0.995, length.out = 500)) |> mutate(fill_value = y)

Figure_5C <- ggplot(ce_depth$DEPTH, aes(x = DEPTH)) +
  geom_raster(data = gradient_data, aes(x = x, y = y, fill = fill_value), interpolate = TRUE) +
  scale_fill_gradient(low = "#FF6666", high = "#9999FF") +
  geom_ribbon(data = poly_low, aes(ymin = ymin, ymax = y), fill = "white") +
  geom_ribbon(data = poly_upp, aes(ymin = y, ymax = ymax), fill = "white") +
  geom_line(aes(y = gam_est), size = 1, color = "black") +
  geom_line(aes(y = gam_low), size = 0.5, color = "black") +
  geom_line(aes(y = gam_upp), size = 0.5, color = "black") +
  geom_point(data = merged_sf, aes(y = Mf_bin, x = DEPTH, fill = Mf_bin), shape = 21, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + theme_classic() +
  labs(x = "Depth", y = "") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text.x     = element_text(size = 16),
        axis.text.y     = element_blank(),      
        axis.ticks.y    = element_blank(),    
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")  

## Conditional effect depth
ce_chl <- conditional_effects(model, effects = "chl") 
ce_chl$chl$gam_low = predict(mgcv::gam(lower__ ~ s(chl, k = 10), data = ce_chl$chl, family = gaussian()), 
                                       newdata = ce_chl$chl)
ce_chl$chl$gam_upp = predict(mgcv::gam(upper__ ~ s(chl, k = 10), data = ce_chl$chl, family = gaussian()), 
                                       newdata = ce_chl$chl)
ce_chl$chl$gam_est = predict(mgcv::gam(estimate__ ~ s(chl, k = 10), data = ce_chl$chl, family = gaussian()), 
                                       newdata = ce_chl$chl)
poly_low <- ce_chl$chl %>% select(chl, gam_low) %>% rename(y = gam_low) %>% arrange(chl) %>%
  mutate(ymin = 0) %>% select(chl, y, ymin)
poly_upp <- ce_chl$chl %>% select(chl, gam_upp) %>% rename(y = gam_upp) %>% arrange(chl) %>%
  mutate(ymax = 1) %>% select(chl, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_chl$chl$chl)+0.001, max(ce_chl$chl$chl)-0.001, length.out = 500), 
                             y = seq(0.01, 0.995, length.out = 500)) |> mutate(fill_value = y)

Figure_5D <- ggplot(ce_chl$chl, aes(x = chl)) +
  geom_raster(data = gradient_data, aes(x = x, y = y, fill = fill_value), interpolate = TRUE) +
  scale_fill_gradient(low = "#FF6666", high = "#9999FF") +
  geom_ribbon(data = poly_low, aes(ymin = ymin, ymax = y), fill = "white") +
  geom_ribbon(data = poly_upp, aes(ymin = y, ymax = ymax), fill = "white") +
  geom_line(aes(y = gam_est), size = 1, color = "black") +
  geom_line(aes(y = gam_low), size = 0.5, color = "black") +
  geom_line(aes(y = gam_upp), size = 0.5, color = "black") +
  geom_point(data = merged_sf, aes(y = Mf_bin, x = chl, fill = Mf_bin), shape = 21, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + theme_classic() +
  labs(x = "Chlorophyll-a concentration", y = "") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text.x     = element_text(size = 16),
        axis.text.y     = element_blank(),      
        axis.ticks.y    = element_blank(),    
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")  

## Conditional effect Functional Diversity
ce_fdiv <- conditional_effects(model, effects = "fdiv_change")
ce_fdiv$fdiv_change$gam_low = predict(mgcv::gam(lower__ ~ s(fdiv_change, k = 10), data = ce_fdiv$fdiv_change, family = gaussian()), 
                                      newdata = ce_fdiv$fdiv_change)
ce_fdiv$fdiv_change$gam_upp = predict(mgcv::gam(upper__ ~ s(fdiv_change, k = 10), data = ce_fdiv$fdiv_change, family = gaussian()), 
                                      newdata = ce_fdiv$fdiv_change)
ce_fdiv$fdiv_change$gam_est = predict(mgcv::gam(estimate__ ~ s(fdiv_change, k = 10), data = ce_fdiv$fdiv_change, family = gaussian()), 
                                      newdata = ce_fdiv$fdiv_change)
poly_low <- ce_fdiv$fdiv_change %>% select(fdiv_change, gam_low) %>% rename(y = gam_low) %>% arrange(fdiv_change) %>%
  mutate(ymin = 0) %>% select(fdiv_change, y, ymin)
poly_upp <- ce_fdiv$fdiv_change %>% select(fdiv_change, gam_upp) %>% rename(y = gam_upp) %>% arrange(fdiv_change) %>%
  mutate(ymax = 1) %>% select(fdiv_change, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_fdiv$fdiv_change$fdiv_change), max(ce_fdiv$fdiv_change$fdiv_change), length.out = 500), 
                             y = seq(0.01, 0.99, length.out = 500)) |> mutate(fill_value = y)

Figure_5E <- ggplot(ce_fdiv$fdiv_change, aes(x = fdiv_change)) +
  geom_raster(data = gradient_data, aes(x = x, y = y, fill = fill_value), interpolate = TRUE) +
  scale_fill_gradient(low = "#FF6666", high = "#9999FF") +
  geom_ribbon(data = poly_low, aes(ymin = ymin, ymax = y), fill = "white") +
  geom_ribbon(data = poly_upp, aes(ymin = y, ymax = ymax), fill = "white") +
  geom_line(aes(y = gam_est), size = 1, color = "black") +
  geom_line(aes(y = gam_low), size = 0.5, color = "black") +
  geom_line(aes(y = gam_upp), size = 0.5, color = "black") +
  geom_point(data = merged_sf, aes(y = Mf_bin, x = fdiv_change, fill = Mf_bin), shape = 21, size = 3) +
  scale_x_continuous(expand = c(0, 0)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  annotate("text", x = min(ce_fdiv$FPI_combined$FPI_combined)+0.01, y = 0.52, 
           label = "no change", hjust = 0, vjust = 0, size = 4, fontface = "italic") +
  coord_cartesian(ylim = c(0, 1)) + theme_classic() +
  labs(x = "Functional diversity change", y = "") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text.x     = element_text(size = 16),
        axis.text.y     = element_blank(),      
        axis.ticks.y    = element_blank(),    
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")  

Figure_5 <- Figure_5A + Figure_5B + Figure_5C + Figure_5D + Figure_5E + plot_layout(ncol = 5)

#### Export the data  ----
## Figures
ggsave(Figure_5, filename = "Figure_5.png", path = "Outputs/", device = "png", width = 18, height = 4, dpi = 300)  
ggsave(Figure_S5, filename = "Figure_S5.png", path = "Outputs/", device = "png", width = 8, height = 5, dpi = 300)  