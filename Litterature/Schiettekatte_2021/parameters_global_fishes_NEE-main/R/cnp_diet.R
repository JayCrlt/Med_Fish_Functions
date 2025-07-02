load_troph <- function(){ ## code to prepare `troph` dataset
      readr::read_csv(here::here("data", "extrapolation_trophic_guilds.csv")) %>%
      dplyr::mutate(genspe = as.character(species)) %>%
      dplyr::select(-species) %>%
      dplyr::mutate(genspe = case_when(genspe == "Chlorurus_sordidus" ~ "Chlorurus_spilururs",
                            TRUE~ genspe))

}

wrangle_dietcnp <- function(data, troph) {
  data %>%
   dplyr::mutate(genspe = gsub(" ", "_", species)) %>%
   dplyr::select(-species, -family) %>%
   left_join(troph) %>%
   # add aulostomus
   mutate(trophic_guild_predicted = case_when(genspe == "Aulostomus_chinensis" ~ (4),
                                              TRUE~ trophic_guild_predicted)) %>%
   mutate(family = case_when(genspe == "Aulostomus_chinensis" ~ "Aulostomidae",
                             TRUE~ as.character(family))) %>%
   filter(key == "AE1", !is.na(n), !is.na(p), !is.na(family), !is.na(trophic_guild_predicted))  %>%
   # remove ouliers per trophic group
   group_by(trophic_guild_predicted) %>%
   filter(abs((p - mean(p))/sd(p)) < 3 ) %>%
   filter(abs((n - mean(n))/sd(n)) < 3 ) %>%
   filter(abs((c - mean(c))/sd(c)) < 3 ) %>%
   dplyr::select(family, genspe,
                 trophic_guild_predicted,
                 id, c, n, p, location, sitename, lat, long,
                 habitat, date, time, weight, tl, sl, sex, year) %>%
   group_by(genspe, location, trophic_guild_predicted) %>%
   mutate(freq = n()) %>%
   ungroup() %>%
   filter(freq>3) %>%
   select(-freq) %>%
   mutate(trophic_guild_predicted = as_factor(trophic_guild_predicted))
}



#' Run model to predict c, n, and p content per trophic guild
#'
#' @param data
#' @param control
#' @param ...
#'
#' @return
#' @export
#' @import brms
#'
#' @examples
#' # run_dietmod()
run_dietmod <- function(data, troph,
                        control = list(adapt_delta = 0.9),
                        ...){
  fit <-  brms::brm(
    mvbind(c, n, p) ~   0 + trophic_guild_predicted + (1|genspe),
    data = data,
    backend= "cmdstanr", threads = threading(8),
    ...
  )

  fit
}

extrapolate_dietcnp <- function(cnpmod, troph, species){
  cnpdiet_fit <- fitted(cnpmod,
                        newdata = data.frame(trophic_guild_predicted =
                                               as.factor(c(1,2,3,4,5,6,7,8))),
                        re_formula =
                          "mvbind(c, n, p) ~   0 + trophic_guild_predicted",
                        ndraws = 1000, summary = FALSE)

  cnpdietfit_df <- lapply(1:8, function(x){
    cnpdiet_fit[,x,] %>% as.data.frame() %>%
      mutate(trophic_guild_predicted = x, draw = 1:1000)}) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(p>0)

  result <- troph %>%
    dplyr::filter(genspe %in% species) %>%
    select(family, genspe, intersect(ends_with("_m"), starts_with("p"))) %>%
    tidyr::pivot_longer(cols = intersect(ends_with("_m"), starts_with("p")),
                 names_to = "pred_cat", values_to = "probability") %>%
    dplyr::mutate(pred_cat = as.factor(pred_cat)) %>%
    dplyr::mutate(pred_cat = fct_recode(pred_cat,
                                 "1" = "p1_m",
                                 "2" = "p2_m",
                                 "3" = "p3_m",
                                 "4" = "p4_m",
                                 "5" = "p5_m",
                                 "6" = "p6_m",
                                 "7" = "p7_m",
                                 "8" = "p8_m")) %>%
    dplyr::right_join(mutate(cnpdietfit_df,
                      pred_cat =
                        as.factor(trophic_guild_predicted))) %>%
    dplyr::group_by(family, genspe, draw) %>%
    dplyr::summarise(c = sum(probability * c)/100,
              n = sum(probability * n)/100,
              p = sum(probability * p)/100) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(family, genspe) %>%
    dplyr::summarise(c_m = mean(c), c_sd = sd(c),
              n_m = mean(n), n_sd = sd(n),
              p_m = mean(p), p_sd = sd(p)) %>%
    dplyr::left_join(select(troph, genspe, trophic_guild_predicted)) %>%
    dplyr::select(family, genspe,
           Dc_m = c_m, Dc_sd = c_m,
           Dn_m = n_m, Dn_sd = n_m,
           Dp_m = p_m, Dp_sd = p_m)
  result
}
