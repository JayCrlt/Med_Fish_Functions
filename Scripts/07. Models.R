#### Setting up          ----
rm(list = ls()) ; options(warn = -1)
library("rfishbase") ; library("phylosem") ; library("tidyverse") ; library('readxl') ; library("scales")
library("fishtree") ; library("geiger") ; library("ape") ; library("Rphylopars") ; library("brms") ; library("sf")
library("ggridges") ; library("patchwork") ; library("fishflux") ; library("leaflet") ; library("leaflet.extras")
library("rnaturalearth") ; library("rnaturalearthdata") ; library("broom") ; library("gstat") ; library("forcats")

## Download data
Env <- read.delim("Data/HAULS_DB_ENV_HINDCAST_COMPLETE_2024.csv", sep = ";")

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
         biom_bin = ifelse(Biomass_slope > 0, 1, 0))

model = brms::brm(formula = Mf_bin ~ botTemp_anom + DEPTH + chl + FPI_combined + (1|GSA), 
                  data = merged_sf, cores = 3, chains = 3, iter = 10000,
                  family = bernoulli(link = "logit"))
(post_summary <- as.data.frame(fixef(model)))
bayes_R2(model)

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
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
    annotate("text", x = min(ce_FPI$FPI_combined$FPI_combined)+0.01, y = 0.52, 
             label = "no change", hjust = 0, vjust = 0, size = 4, fontface = "italic") +
    coord_cartesian(ylim = c(0, 1)) + theme_classic() +
    labs(x = "Fish Pressure Index", y = "Multifunctional Index change") +
     theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
           plot.title      = element_text(size = 20),
           axis.title      = element_text(size = 18),
           axis.text.x     = element_text(size = 16),
           axis.text.y     = element_blank(),      
           axis.ticks.y    = element_blank(),    
           legend.title    = element_text(size = 18),
           legend.text     = element_text(size = 16),
           legend.position = "none")  
  
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
  labs(x = "Bottom Temperature Anomaly", y = "") +
  theme(panel.border    = element_rect(color = "black", fill = NA, size = 1),
        plot.title      = element_text(size = 20),
        axis.title      = element_text(size = 18),
        axis.text.x     = element_text(size = 16),
        axis.text.y     = element_blank(),      
        axis.ticks.y    = element_blank(),    
        legend.title    = element_text(size = 18),
        legend.text     = element_text(size = 16),
        legend.position = "none")  

ce_sbta <- conditional_effects(model, effects = "DEPTH") 
ce_sbta$DEPTH$gam_low = predict(mgcv::gam(lower__ ~ s(DEPTH, k = 10), data = ce_sbta$DEPTH, family = gaussian()), 
                                       newdata = ce_sbta$DEPTH)
ce_sbta$DEPTH$gam_upp = predict(mgcv::gam(upper__ ~ s(DEPTH, k = 10), data = ce_sbta$DEPTH, family = gaussian()), 
                                       newdata = ce_sbta$DEPTH)
ce_sbta$DEPTH$gam_est = predict(mgcv::gam(estimate__ ~ s(DEPTH, k = 10), data = ce_sbta$DEPTH, family = gaussian()), 
                                       newdata = ce_sbta$DEPTH)
poly_low <- ce_sbta$DEPTH %>% select(DEPTH, gam_low) %>% rename(y = gam_low) %>% arrange(DEPTH) %>%
  mutate(ymin = 0) %>% select(DEPTH, y, ymin)
poly_upp <- ce_sbta$DEPTH %>% select(DEPTH, gam_upp) %>% rename(y = gam_upp) %>% arrange(DEPTH) %>%
  mutate(ymax = 1) %>% select(DEPTH, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_sbta$DEPTH$DEPTH)+0.001, max(ce_sbta$DEPTH$DEPTH)-0.001, length.out = 500), 
                             y = seq(0.01, 0.995, length.out = 500)) |> mutate(fill_value = y)

Figure_5C <- ggplot(ce_sbta$DEPTH, aes(x = DEPTH)) +
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

ce_sbta <- conditional_effects(model, effects = "chl") 
ce_sbta$chl$gam_low = predict(mgcv::gam(lower__ ~ s(chl, k = 10), data = ce_sbta$chl, family = gaussian()), 
                                       newdata = ce_sbta$chl)
ce_sbta$chl$gam_upp = predict(mgcv::gam(upper__ ~ s(chl, k = 10), data = ce_sbta$chl, family = gaussian()), 
                                       newdata = ce_sbta$chl)
ce_sbta$chl$gam_est = predict(mgcv::gam(estimate__ ~ s(chl, k = 10), data = ce_sbta$chl, family = gaussian()), 
                                       newdata = ce_sbta$chl)
poly_low <- ce_sbta$chl %>% select(chl, gam_low) %>% rename(y = gam_low) %>% arrange(chl) %>%
  mutate(ymin = 0) %>% select(chl, y, ymin)
poly_upp <- ce_sbta$chl %>% select(chl, gam_upp) %>% rename(y = gam_upp) %>% arrange(chl) %>%
  mutate(ymax = 1) %>% select(chl, y, ymax)
gradient_data <- expand.grid(x = seq(min(ce_sbta$chl$chl)+0.001, max(ce_sbta$chl$chl)-0.001, length.out = 500), 
                             y = seq(0.01, 0.995, length.out = 500)) |> mutate(fill_value = y)

Figure_5D <- ggplot(ce_sbta$chl, aes(x = chl)) +
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

Figure_5 <- Figure_5A + Figure_5B + Figure_5C + Figure_5D + plot_layout(ncol = 4)

#### Export the data  ----
## Figures
ggsave(Figure_5, filename = "Figure_5.png", path = "Outputs/", device = "png", width = 15, height = 4, dpi = 300)  