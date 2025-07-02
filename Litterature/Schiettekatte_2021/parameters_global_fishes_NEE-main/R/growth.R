run_kmax_model <- function(fishtree, oto) {
  # Prerequesite for phylo brms is the covariance matrix (following the example in the brms phylo vignette)
  # --> Do this for each tree
  corphy_list <- parallel::mclapply(fishtree, function(x){
    phy2 <- x %>% ape::chronoMPL()
    inv.phylo <- MCMCglmm::inverseA(phy2, nodes = "TIPS", scale = TRUE)
    A <- solve(inv.phylo$Ainv)
    rownames(A) <- rownames(inv.phylo$Ainv)
    A[sort(rownames(A)), sort(rownames(A))]
    return(A)
  }, mc.cores = 40)

  # Transform list of covariance matrices to 3d array
  corphy_array <- simplify2array(corphy_list)

  # Summarize these matrices
  corphy_m <- apply(corphy_array, 1:2, mean)
  corphy_sd <- apply(corphy_array, 1:2, sd)

  ### subset only for oto data
  tt <- which(rownames(corphy_m) %in% unique(oto$species))
  corphy_oto <- corphy_m[tt, tt]
  corphy_oto <- corphy_oto[sort(rownames(corphy_oto)), sort(rownames(corphy_oto))]

  oto <- filter(oto, species %in% rownames(corphy_m)) %>%
    filter(method == "Otolith rings")

  oto$phylo <- as.factor(oto$species)

  ##### model  #####

  fit_kmax <- brm(
    log(kmax) ~ log(sizemax) + sstmean + (1|phylo), data = oto,
    family = gaussian(), cov_ranef = list(phylo = corphy_oto), cores = 4,
    backend = "cmdstanr", threads = threading(8)
  )

  pp_check(fit_kmax)
  bayes_R2(fit_kmax)

  # phylogenetic signal
  hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
  (hyp <- hypothesis(fit_kmax, hyp, class = NULL))
  hyp

  return(fit_kmax)

}


##### extrapolate #####
extrapolate_kmax <- function(fit_kmax, fishtree, sp_sst){

  r_phylo <- tidybayes::spread_draws(fit_kmax, r_phylo[species,Intercept])

  b <- tidybayes::spread_draws(fit_kmax, b_Intercept, b_logsizemax, b_sstmean)

  # each draw
  # round sst to be reasonable
  sp_sst$sst <- round(sp_sst$sst)
  sp_sst <- unique(sp_sst)

  all_draw_kmax <- lapply(1:1000, function(x){

  print(paste("draw", x))
  draw = x

  b_sub <- filter(b, .draw == draw)
  b_int <- b_sub$b_Intercept
  b_logsize <- b_sub$b_logsizemax
  b_sst <- b_sub$b_sstmean

  r_spec <- dplyr::filter(r_phylo, .draw == draw)
  st1 <- r_spec$r_phylo
  names(st1) <- r_spec$species

  rphy_list <- mclapply(1:100, function(x){
    r_phy <- picante::phyEstimate(fishtree[[x]], st1)
    return(r_phy)
  }, mc.cores = 40)

  rphy_l <- lapply(rphy_list, function(x){ return(x[,1])})

  tphy_array <- simplify2array(rphy_l)
  dim(tphy_array)

  rphy_m <- data.frame(
    species = rownames(rphy_list[[1]]),
    r_phylo = apply(tphy_array, 1, mean, na.rm = TRUE))
  rphy_1 <- data.frame(
    species = names(st1),
    r_phylo = st1
  )

  kmax_ex <- rbind(rphy_m, rphy_1) %>% group_by(species) %>% summarise_all(mean) %>%
    filter(species %in% fishtree[[1]]$tip.label) %>% right_join(sp_sst) %>% unique() %>%
    mutate(b_logsize = b_logsize, b_sst = b_sst, b_int = b_int) %>%
    mutate(logkmax = r_phylo + b_logsize * log(max_size) + b_sst * sst + b_int) %>% mutate(kmax = exp(logkmax), draw = draw)

  return(kmax_ex)

  })

  kmaxtall <-plyr::ldply(all_draw_kmax)

  kmax_predict <- plyr::ldply(all_draw_kmax) %>% group_by(Family, species, max_size, sst) %>%
    summarise(kmax_m = mean(kmax), kmax_sd = sd(kmax), r_phylo_m = mean(r_phylo), r_phylo_sd = sd(r_phylo))

  return(kmax_predict)

}
