library(ape)
library(fishtree)
library(brms)
library(purrr)
library(tidybayes)
library(rstan)
library(parallel)
library(Hmisc)
library(Rphylopars)
library(tidyverse)
library(fishtree)
library(Rphylopars)
library(Hmisc)
library(ggridges)
library(ggstance)
library(rfishbase)
library(magrittr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(drake)

options(mc.cores = parallel::detectCores())



###### load global fish data ########

load(here::here("data","fishtree_glob.RData")) # source: fishtree
