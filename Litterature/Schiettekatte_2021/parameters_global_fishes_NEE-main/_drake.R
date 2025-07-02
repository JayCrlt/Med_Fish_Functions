source("R/packages.R")  # loads packages
source("R/cnp_body.R")
source("R/cnp_diet.R")
source("R/growth.R")
source("R/metabolism.R")
source("R/plan.R")  


future::plan(future.callr::callr)

config <- drake_config(plan,
                       lock_envir = FALSE,
                       garbage_collection = TRUE,
                       memory_strategy = "autoclean",
                       history = FALSE)

