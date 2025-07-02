
wrangle_resp <- function(fdata){
  tax <- rfishbase::load_taxa() %>%
    dplyr::select(Species, Family) %>%
    as.data.frame()
  fdata_sp <- as.data.frame(table(fdata$Species))
  fdata_sp <- filter(fdata_sp, Freq > 2)
  colnames(fdata_sp) <- c("Species", "freq")
  fdata_sp <- left_join(fdata_sp, tax)


  fdata <- filter(fdata, Species %in% fdata_sp$Species)%>%
    inner_join(tax) %>% tidyr::drop_na()

  fdata$mass <- fdata$Weight..kg.*1000
  fdata$SMR_gd <- fdata$SMR * (12/32) * 0.001 * 24
  fdata$MMR_gd <- pmax(fdata$FirstMR,fdata$MaxMR) * (12/32) * 0.001 * 24


  fdata_long <- tidyr::gather(fdata,"mr_type", "mr", SMR_gd, MMR_gd)
  fdata_long$mr_type <- as.factor(fdata_long$mr_type)
  fdata_long$Species <- as.factor(fdata_long$Species)
  fdata_long$Family <- as.factor(fdata_long$Family)
  fdata_long <- filter(fdata_long, mr>0)

  fdata_long
}

run_model_resp <- function(respdata){

  fit_mr <- brm(log10(mr) ~ 1 + log10(mass) + (0 + log10(mass)|Family) + (0 + log10(mass)|Species) +
                  (1|Family:mr_type) + (1|Species:mr_type) + (0 + log10(mass)|Family:mr_type) +
                  (0 + log10(mass)|Species:mr_type), data = respdata,
                control = list(adapt_delta = 0.9),#, max_treedepth = 12),
                prior = set_prior(class = "b", "log10mass", prior = "normal(0.8, 0.5)"), cores = 4,
                warmup = 2000, iter = 4000, backend= "cmdstanr", threads = threading(8))
  pp_check(fit_mr)

  return(fit_mr)

}


exp10 <- function(x) {
  10^x
}

predict_resp_family <- function(fit_mr) {
  draw_f<- fit_mr %>%
    spread_draws(r_Family[Family,par]) %>% as.data.frame()

  draw_fmr<- fit_mr %>%
    spread_draws(`r_Family:mr_type`[Family,par]) %>% as.data.frame() %>%
    separate(col = "Family", into = c("Family", "mr", "unit"), sep = "_")

  draw_b <-  fit_mr %>%
    spread_draws(b_log10mass) %>% as.data.frame()

  draw_int <-  fit_mr %>%
    spread_draws(b_Intercept) %>% as.data.frame()

  draws_b <- left_join(draw_f, draw_b) %>% left_join(filter(draw_fmr, par == "log10mass")) %>%
    mutate(b = b_log10mass + r_Family + `r_Family:mr_type`)
  draws_int <- left_join(filter(draw_fmr, par == "Intercept"), draw_int)  %>%
    mutate(int = b_Intercept + `r_Family:mr_type`)


  alpha <- draws_b %>% group_by(Family, mr) %>% summarise(alpha = mean(b), alpha_sd = sd(b))

  B0 <- draws_int %>% group_by(Family, mr) %>% summarise(B0 = mean(exp10(int)), B0_sd = sd(exp10(int)))
  B0_smr <- filter(B0, mr == "SMR")
  B0_mmr <- filter(B0, mr == "MMR")
  far <- (3 * B0_smr$B0 + B0_mmr$B0) /(4 * B0_smr$B0)
  far <- data.frame(Family = unique(B0$Family), theta = far)
  B0_mmr$B0/B0_smr$B0

  metpar <- left_join(alpha, B0) %>% left_join(far)
  metpar_smr <- filter(metpar, mr == "SMR") %>% select(-mr)

  metpar$B0 * exp(0.59 / 8.62e-5 * ( (1 / (28 + 273.15)) - (1 / (38 + 273.15))))

  return(metpar_smr)


}
