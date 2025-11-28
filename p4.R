# PART 4 (resolve) ------------------------------------------------------------------

# STEP 1: Resolve trends & occupancy for all selected species
# Run:
# - after above steps (occupancy)
# Requires:
# - tidyverse, tictoc, sf, VGAM, writexl
# - data files:
#   - fullspecieslist.csv
#   - trends/trendsX.csv for whole country and individual mask versions
# Outputs: several

library(tidyverse)
library(glue)
library(tictoc)
# for parallel iterations
library(furrr)
library(parallel)

interannual_update = FALSE
load("00_data/analyses_metadata.RData")
cur_mask <- "none"

species_to_process <- c("Black-crested Bulbul",
                        "Malabar Whistling-Thrush",
                        "Greater Racket-tailed Drongo",
                        "Streaked Spiderhunter",
                        "Gray Francolin",
                        "Spotted Dove",
                        "Greater Coucal",
                        "Red Junglefowl",
                        "Little Grebe",
                        "Little Cormorant",
                        "Montagu's Harrier",
                        "Common Chiffchaff",
                        "Tawny Pipit",
                        "Nilgiri Wood-Pigeon",
                        "Black-throated Thrush",
                        "Northern Lapwing",
                        "Little Swift",
                        "Black-winged Kite")

# cur_metadata <- get_metadata(cur_mask)
# 
# # read paths
# base_path <- cur_metadata$FULLSPECLIST.PATH
# 
# fullspecieslist <- read.csv(base_path) %>%
# rename(proprange25km.latestyear = proprange25km2022)
# 
# write.csv(fullspecieslist, "01_analyses_full/fullspecieslist.csv")
  
tic.clearlog()
tic("Resolved trends & occupancy for full country")
source("00_scripts/resolve_trends_and_occupancy.R")
toc(log = TRUE, quiet = TRUE) 
tic.log()

# STEP 2: Classify using trends and range status, and generate necessary outputs
# Run:
# - after above steps (P4, S1)
# Requires:
# - tidyverse, tictoc, writexl
# - data files:
#   - specieslists.RData
#   - trends/trendsX.csv for whole country and individual mask versions
#   - X/SoIB_main_wocats.csv
# Outputs: several

tic.clearlog()
tic("Finished classifying and summarising for full country") # 2 min
source("00_scripts/classify_and_summarise.R")
toc(log = TRUE)
tic.log()    
