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
    dplyr::select(-Estimated, -all_of(trait_col)) |> rename(!!trait_name := paste0(trait_name, "_final"))
  return(result)}


# Imputational function FLUXGLOB
impute_trait_phylopars_FG <- function(data, trait_col, phy) {
  data_phy <- data |> dplyr::mutate(Species_phy = gsub(" ", "_", Species), 
                                     Trait_orig = .data[[trait_col]]) |> dplyr::ungroup()
  trait_data <- data_phy |> dplyr::filter(!is.na(Trait_orig)) |>
    dplyr::transmute(Species = Species_phy, Trait = as.numeric(Trait_orig)) |>
    dplyr::distinct(Species, .keep_all = TRUE)
  matched_species <- intersect(trait_data$Species, phy$tip.label)
  trait_data_clean <- trait_data |> dplyr::filter(Species %in% matched_species) |>
    dplyr::transmute(species = Species, Trait = Trait) |> as.data.frame()
  if (nrow(trait_data_clean) < 2) {stop(paste0("Not enough data for trait: ", trait_col))}
  fit <- Rphylopars::phylopars(trait_data = trait_data_clean, tree = phy)
  imputed_matrix <- fit$anc_recon[phy$tip.label, , drop = FALSE]
  estimates <- tibble::tibble(Species_phy = rownames(imputed_matrix), Estimated = imputed_matrix[, "Trait"])
  trait_name <- gsub("[^[:alnum:]]", "", trait_col)
  result <- data_phy |>  dplyr::left_join(estimates, by = "Species_phy") |>
    dplyr::mutate(!!trait_col := dplyr::coalesce(Trait_orig, Estimated),
                  !!paste0(trait_name, "_type") := dplyr::case_when(
                    !is.na(Trait_orig) ~ "measured",
                    is.na(Trait_orig) & !is.na(Estimated) ~ "phylogenetic", TRUE ~ "missing")) |>
    dplyr::select(-Species_phy, -Estimated, -Trait_orig)
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
  df <- df |>
    group_by(Genus) |>
    mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE),!!trait_sym)) |>
    ungroup() |> group_by(Family) |>
    mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE),!!trait_sym)) |>
    ungroup() |> mutate(!!trait_sym := ifelse(is.na(!!trait_sym), mean(!!trait_sym, na.rm = TRUE), !!trait_sym))
  return(df)}


# Function to fill NAs hierarchically for one trait FLUXGLOB
fill_trait_hierarchy_FG <- function(df, trait) {
  fill_level <- function(data, level) {
  data |> group_by(across(all_of(level))) |> mutate("{trait}" := ifelse(is.na(.data[[trait]]) & 
                                                              !is.nan(mean(.data[[trait]], na.rm = TRUE)),
          mean(.data[[trait]], na.rm = T), .data[[trait]])) |> ungroup()}
  df |> fill_level("Genus") |> fill_level("Subfamily") |> fill_level("Family") |> 
    mutate( "{trait}" := ifelse(is.na(.data[[trait]]), mean(.data[[trait]], na.rm = TRUE), .data[[trait]]))}


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


# Temperature FLUXGLOB
impute_temp_FG <- function(month, year, xy) {
  cand <- temp_available[temp_available$Month == month, ]
  if (nrow(cand) == 0) cand <- temp_available
  cand_coords <- cbind(as.numeric(cand$X), as.numeric(cand$Y))
  nn <- FNN::get.knnx(cand_coords, matrix(xy, nrow = 1), k = min(20, nrow(cand)))
  neigh <- cand[nn$nn.index[1, ], ]
  neigh$year_diff <- abs(neigh$Year - year)
  neigh |> arrange(year_diff) |> slice(1) |> pull(SBTemp)}


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


# Clean numeric
clean_numeric <- function(x) as.numeric(gsub(",", ".", x))


# label Extended Figure
make_label <- function(id, eco) {
  map2_chr(id, eco, function(i, e) {
    eco_wrapped <- str_wrap(e, width = 25)
    n_lines <- str_count(eco_wrapped, "\n") + 1
    eco_html <- gsub("\n", "<br>", eco_wrapped)
    if (n_lines == 1) { paste0("<b>", i, ".</b> ", eco_html, "<br>&nbsp;") } else {
      paste0("<b>", i, ".</b> ", eco_html)}})}