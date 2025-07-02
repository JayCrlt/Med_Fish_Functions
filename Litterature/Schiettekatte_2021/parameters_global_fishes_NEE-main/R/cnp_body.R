run_cnp_model <- function(cnptree, cnpdata) {

  # covariance matric based on phylogeny
  corphy_list <- lapply(cnptree, function(x){
  phy <-ape::chronoMPL(x)
  inv.phylo <- MCMCglmm::inverseA(phy, nodes = "TIPS", scale = TRUE)
  A <- solve(inv.phylo$Ainv)
  rownames(A) <- rownames(inv.phylo$Ainv)
  A <- A[sort(rownames(A)), sort(rownames(A))]
  return(A)
})

corphy_array <- simplify2array(corphy_list)
corphy_m <- apply(corphy_array, 1:2, mean)
cnpdata$phylo <- cnpdata$species

# run model
fit_cnp <- brm(
  mvbind(C, N, P) ~  (1|p|phylo) + (1|species),
  data = cnpdata, family = gaussian(),
  cov_ranef = list(phylo = corphy_m),
  sample_prior = TRUE, chains = 2, cores = 4,
  iter = 2000, warmup = 1000, backend= "cmdstanr", threads = threading(8)
)

hyp <- "sd_phylo__C_Intercept^2 / (sd_phylo__C_Intercept^2 + sigma_C^2) = 0"
(hyp <- hypothesis(fit_cnp, hyp, class = NULL)) # 0.35 (0.24, 0.51)
hyp <- "sd_phylo__N_Intercept^2 / (sd_phylo__N_Intercept^2 +sigma_N^2) = 0"
(hyp <- hypothesis(fit_cnp, hyp, class = NULL)) # 0.53  ( 0.4 ,  0.66)
hyp <- "sd_phylo__P_Intercept^2 / (sd_phylo__P_Intercept^2 +  sigma_P^2) = 0"
(hyp <- hypothesis(fit_cnp, hyp, class = NULL)) # 0.56 (0.42, 0.68)

brms::bayes_R2(fit_cnp)

pp_check(fit_cnp, resp = "C")
pp_check(fit_cnp, resp = "N")
pp_check(fit_cnp, resp = "P")


return(fit_cnp)

}

extract_draws_cnp_phylo <- function(fit, cnpdata){
  cnpdata$phylo <- cnpdata$species
  fitted(fit, newdata = unique(select(cnpdata, species, phylo)),
                      summary = FALSE)
}

extrapolate_bodycnp <- function(cnpdata, fit, spglob, tree, cnp_phylo) {

  test <-
    lapply(1:1000, function(x){
      print(x)
      trait <- data.frame(
        species = unique(cnpdata$species),
        C = cnp_phylo[x, ,1],
        N = cnp_phylo[x, ,2],
        P = cnp_phylo[x, ,3]
      ) %>% unique() %>% group_by(species) %>% summarise_all(mean)
      trait <- left_join(spglob, trait) %>% select(species, C, N, P) %>%
        as.data.frame()
      all <- parallel::mclapply(1:100, function(y){
        phylopars(trait, fishtree[[y]])$anc_recon
        }, mc.cores = 50)
      alls <- simplify2array(all) %>% apply( 1:2, mean)
      return(alls)
    })

  mm <- lapply(test, max) %>% simplify2array()
  torm <-which((mm)>90)


  testa <-  simplify2array(test)
  testa[1,,]
  dim(testa)
  testa <- testa[,,-torm]

  trait_m <- apply(testa, 1:2, median) %>% as.data.frame()
  trait_m$species <- rownames(trait_m)
  trait_m <- trait_m %>% select(species, C, N, P)
  trait_sd <- apply(testa, 1:2, sd) %>% as.data.frame()
  colnames(trait_sd) <- c("C_sd", "N_sd", "P_sd")

  cnp_all <- cbind(trait_m, trait_sd) %>% inner_join(spglob) %>%
    dplyr::select(species, Qc_m = C, Qc_sd = C_sd,
                  Qn_m = N, Qn_sd = N_sd,
                  Qp_m = P, Qp_sd = P_sd)

  cnp_all

}
