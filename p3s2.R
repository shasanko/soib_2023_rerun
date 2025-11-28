# STEP 2: Run occupancy analyses (presence-based and model)
# Run:
# - after p3s1a.R
# - both analyses run only for full country;
# - this is not run for habmasks at all
# Requires:
# - tidyverse, tictoc, glue, parallel, foreach, doParallel
# - data files:
#   - "dataforanalyses.RData" for whole country and individual states
#   - "specieslists.RData" for whole country and individual states
#   - "00_data/SoIB_mapping_2022.csv"
#   - "00_data/grids_sf_nb.RData"
# Outputs: 
# - csv files in occupancy-presence/ 
# - "occupancy-model/chunk_X.csv" for whole country and individual states

library(tidyverse)
library(glue)
library(tictoc)
# for parallel iterations
library(furrr)
library(parallel)

interannual_update = FALSE
source("00_scripts/00_functions.R")
load("00_data/analyses_metadata.RData")

# species_to_process <- c("Black-crested Bulbul",
#                         "Malabar Whistling-Thrush",
#                         "Greater Racket-tailed Drongo",
#                         "Streaked Spiderhunter",
#                         "Gray Francolin",
#                         "Spotted Dove",
#                         "Greater Coucal",
#                         "Red Junglefowl",
#                         "Little Grebe",
#                         "Little Cormorant",
#                         "Montagu's Harrier",
#                         "Common Chiffchaff",
#                         "Tawny Pipit",
#                         "Nilgiri Wood-Pigeon",
#                         "Black-throated Thrush",
#                         "Northern Lapwing",
#                         "Little Swift",
#                         "Black-winged Kite")

species_to_process <- c("Black-crested Bulbul",
                        "Malabar Whistling-Thrush",
                        "Streaked Spiderhunter",
                        "Gray Francolin"
                        )
 
# full country
cur_mask <- "none"
source("00_scripts/run_species_occupancy-setup.R")
tic("Ran presence-based occupancy")
source("00_scripts/run_species_occupancy-presence.R")
toc()  
tic("Ran modelled occupancy")
source("00_scripts/run_species_occupancy-model.R")
toc()  


# occupancy not run for hab masks at all. both presence and modelled data pulled from full-country.



