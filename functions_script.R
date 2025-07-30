# Imputational function
impute_trait_phylopars <- function(data, trait_col, phy) {
  trait_vector <- data |>
    filter(!is.na(.data[[trait_col]])) |>
    transmute(Species = gsub(" ", "_", Species), Trait = .data[[trait_col]]) |>
    tibble::deframe()
  trait_data <- data.frame(Species = names(trait_vector), Trait = trait_vector)
  trait_data$Species <- gsub(" ", "_", trait_data$Species)
  matched_species <- intersect(trait_data$Species, phy$tip.label)
  trait_data_clean <- trait_data |> filter(Species %in% matched_species)
  fit <- phylopars(trait_data = trait_data_clean |> rename(species = Species), tree = phy)
  imputed_matrix <- fit$anc_recon[phy$tip.label, , drop = FALSE]
  trait_name <- gsub("[^a-zA-Z]", "", trait_col)  
  estimates <- tibble(Species = rownames(imputed_matrix), Estimated = imputed_matrix[, "Trait"])
  result <- data |> 
    left_join(estimates, by = "Species") |>
    mutate(
      !!paste0(trait_name, "_final") := ifelse(is.na(.data[[trait_col]]), Estimated, .data[[trait_col]]),
      !!paste0(trait_name, "_type") := ifelse(is.na(.data[[trait_col]]), "imputed", "measured")) |>
    select(-Estimated, -all_of(trait_col)) |> rename(!!trait_name := paste0(trait_name, "_final"))
  return(result)}


# Sampling
sample_once <- function(df, n_sample) {
  df |> group_by(YEAR) |> slice_sample(n = n_sample) |>
    mutate(HAUL_DURATION = as.numeric(HAUL_DURATION), DISTANCE = as.numeric(DISTANCE), TOTAL_WEIGHT_DIST_KG_KM_H = 
             (TOTAL_WEIGHT_IN_THE_HAUL / 1000) / (DISTANCE / 1000) / (HAUL_DURATION * 60)) |>
    group_by(YEAR) |> summarise(weight_std = sum(TOTAL_WEIGHT_DIST_KG_KM_H, na.rm = T)) |> ungroup()}


# Function to fill NAs hierarchically for one trait
fill_trait_hierarchy <- function(df, trait) {
  trait_sym <- sym(trait)
  df <- df %>%
    group_by(Genus) %>%
    mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE),!!trait_sym)) %>%
    ungroup() %>% group_by(Family) %>%
    mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE),!!trait_sym)) %>%
    ungroup() %>% mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE), !!trait_sym))
  return(df)}