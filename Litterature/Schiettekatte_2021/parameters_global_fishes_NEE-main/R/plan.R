plan <- drake_plan(

  ##### set up output #####
  out_folder = dir.create("output"),

  # global species list
  spglob = readr::read_csv(here::here("data", "species_list_tl.csv")),
  sp_sst = readr::read_csv(here::here("data", "sp_sst.csv")),
  troph = load_troph(),
  fishtree_glob = load(here::here("data","fishtree_glob.RData")), # source: fishtree,


  ##### body cnp #####
  bodycnp = readr::read_csv(here::here("data","cnpdata_combined.csv")),
  fishtree_bodycnp = fishtree_complete_phylogeny(species = unique(bodycnp$species)),
  model_bodycnp = run_cnp_model(fishtree_bodycnp, bodycnp),
  draws_cnp_phylo = extract_draws_cnp_phylo(model_bodycnp, bodycnp),
  bodycnp_extrap = extrapolate_bodycnp(bodycnp, model_bodycnp, spglob, fishtree_glob, draws_cnp_phylo),
  out_bodycnp_extrap = write_csv(bodycnp_extrap, here::here("output", "bodycnp_extrapolated.csv")),

  ##### growth #####
  oto = readr::read_csv(here::here("data","kmax_combined.csv")),
  model_kmax = run_kmax_model(fishtree, oto),
  kmax_extrap = extrapolate_kmax(model_kmax, fishtree, sp_sst),
  out_kmax_extrap = write_csv(kmax_extrap, here::here("output", "kmax_predict.csv")),

  ##### diet cnp #####
  dietcnp_raw = readr::read_csv(here::here("data","cnp_diet.csv")),
  dietcnp = wrangle_dietcnp(dietcnp_raw, troph),
  model_dietcnp = run_dietmod(dietcnp, troph),
  dietcnp_extrap = extrapolate_dietcnp(model_dietcnp, troph, spglob),
  out_dietcnp_extrap = write_csv(dietcnp_extrap, here::here("output", "dietcnp_extrapolated.csv")),

  ##### metabolic rate #####
   data_resp = readr::read_csv(here::here("data", "fishdata_respiro_2019.csv")),
   resp_clean = wrangle_resp(data_resp),
   model_resp = run_model_resp(resp_clean),
   resp_family = predict_resp_family(model_resp),
   out_resp_family = write_csv(resp_family, here::here("output", "metpar_fam_smr.csv"))
)
