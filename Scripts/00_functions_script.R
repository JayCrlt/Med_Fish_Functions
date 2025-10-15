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


# Ranking 1% functions
top_1_percent <- function(x) {
  threshold <- quantile(x, 0.99, na.rm = TRUE)
  x >= threshold}


# Temperature
impute_temp <- function(hex_id, date, month) {
  candidates <- temp_available[temp_available$HEX_ID == hex_id, ]
  if (nrow(candidates) > 0) {
    same_month <- candidates[candidates$Month == month, ]
    if (nrow(same_month) > 0) {
      same_month$diff <- abs(as.numeric(difftime(same_month$Date, date, units="days")))
      return(same_month$TEMP[which.min(same_month$diff)])}
    candidates$diff <- abs(as.numeric(difftime(candidates$Date, date, units="days")))
    return(candidates$TEMP[which.min(candidates$diff)])}
  my_coords <- hex_lookup[hex_lookup$HEX_ID == hex_id, c("X","Y")]
  nn <- get.knnx(hex_lookup[,c("X","Y")], my_coords, k = nrow(hex_lookup))
  ordered_hex <- hex_lookup$HEX_ID[nn$nn.index[1,]]
  for (h in ordered_hex[-1]) {
    candidates <- temp_available[temp_available$HEX_ID == h, ]
    if (nrow(candidates) > 0) {same_month <- candidates[candidates$Month == month, ]
    if (nrow(same_month) > 0) {same_month$diff <- abs(as.numeric(difftime(same_month$Date, date, units="days")))
    return(same_month$TEMP[which.min(same_month$diff)])}
    candidates$diff <- abs(as.numeric(difftime(candidates$Date, date, units="days")))
    return(candidates$TEMP[which.min(candidates$diff)])}}; return(NA) }


# Percentile
percentile_class <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.99, 0.05, 0.95), na.rm = TRUE)
  case_when(
    x <= q[1] ~ "tail_1%",
    x >= q[2] ~ "top_1%",
    x <= q[3] ~ "tail_5%",
    x >= q[4] ~ "top_5%",
    TRUE ~ "mid")}


# Not in
`%notin%` = Negate(`%in%`)