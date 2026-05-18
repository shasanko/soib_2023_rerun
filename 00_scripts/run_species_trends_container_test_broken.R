library(tidyverse)

library(lme4)
library(VGAM)
library(parallel)

container <- FALSE

reproducible_run <- TRUE

hostname <- paste0(Sys.info()["nodename"],"")

# preparing data for specific mask (this is the only part that changes, but automatically)
cur_metadata <- get_metadata(cur_mask, container)
data_prefix = ""
speclist_path <- paste0(data_prefix, cur_metadata$SPECLISTDATA.PATH)
databins_path <- paste0(data_prefix, cur_metadata$DATA.PATH) # for databins

get_free_ram <- function() {
  #               total        used        free      shared  buff/cache   available
  #Mem:        65572748     1544972    62672624        2976     1995520    64027776
  #Swap:        8388604           0     8388604
  #"available" is how much more RAM we can be use without swapping
  ram <- system("free | awk '/Mem:/ {print $7}'", intern = TRUE)
  ram <- as.integer(ram)*1024
  return(ram)
}

# don't run if no species selected
message(paste("Loading:",speclist_path))
load(speclist_path)
to_run <- (1 %in% specieslist$ht) | (1 %in% specieslist$rt) |
  (1 %in% restrictedspecieslist$ht) | (1 %in% restrictedspecieslist$rt)


# singleyear = interannualupdate
singleyear = FALSE
# not using single year modelling approach, since test runs showed that
# single year models produce notably higher estimates than full-year models

# for the full country analysis, runs are split among multiple systems, and use
# separate subsampled datasets. We need to ensure this information exists.
# else, all 1000 runs are on one system.


source('00_scripts/00_functions.R')

databins_path_metadata <- paste0(databins_path,'-metadata')
message(paste("Loading:",databins_path_metadata))
load(paste(databins_path_metadata))

databins_path_data <- paste0(databins_path,'-data_opt')
tic(paste("Loading:",databins_path_data))
load(paste(databins_path_data)) # will create "data"
toc()

lsa = specieslist %>% filter(!is.na(ht) | !is.na(rt))
listofspecies = c(lsa$COMMON.NAME, restrictedspecieslist$COMMON.NAME)
speclen = length(listofspecies)

# creating new directory if it doesn't already exist
if (!dir.exists(cur_metadata$TRENDS.PATHONLY)) {
  dir.create(cur_metadata$TRENDS.PATHONLY, 
             recursive = T)
}

# Load common data
basedir <- dirname(databins_path)
species_names_path <- paste0(basedir, "/species_names.RData")
timegroups_path <- paste0(basedir, "/timegroups.RData")
message("Loading: ", species_names_path)
load(species_names_path)
message("Loading: ", timegroups_path)
load(timegroups_path)

# run_stats_path <- paste0(dirname(databins_path),'/species_run_stats.RData')
# have_run_stats <- FALSE
# if(file.exists(run_stats_path)) {
#   message(paste("Loading:", run_stats_path))
#   load(run_stats_path)
#   have_run_stats <- TRUE
# } else {
#   message("No run stats, can't optimize")
# }

# delete coulumns gridg2 and gridg4
# data$OBSERVER.ID <- NULL

k <- 1

species_to_process <- c("Gray Francolin")

message("========================================")
message(paste("Starting assignment:", k))
message("========================================")

trends_species_dir <- cur_metadata %>%
  dplyr::summarise(TRENDS.PATH = glue("{TRENDS.PATHONLY}/species_step_through_{k}")) %>%
  as.character()

trends_stats_dir <- paste0(trends_species_dir,"/stats")
# creating new directory if it doesn't already exist
if (!dir.exists(trends_stats_dir)) {
  dir.create(trends_stats_dir,
             recursive = T)
}

# file names for individual files

write_path <- cur_metadata %>%
  dplyr::summarise(TRENDS.PATH = glue("{TRENDS.PATHONLY}trends_step_through_{k}.csv")) %>%
  as.character()

# data_path = cur_metadata %>%
#   dplyr::summarise(SIMDATA.PATH = glue("{data_prefix}{SIMDATA.PATHONLY}data{k}.RData_opt")) %>%
#   as.character()

tictoc::tic(glue("Species trends for {cur_mask}: {k}/{max(cur_assignment)}"))

# read data files for this step
rgid_path <- paste0(dirname(databins_path_data),"/rgids-", k, ".RData")
message(paste("Loading", rgid_path))
load(rgid_path) # loads randomgroupids

# Subset data to match assignment
data_filt <- data[data$group.id %in% randomgroupids, ]

#save(data_filt, file = "01_analyses_full/data_filt_k_1.RData")
#write.csv(data_filt, file = "data_filt_k_1.csv")

# map timegroups to strings
data_filt$timegroups <- timegroups_names$timegroups[data_filt$timegroups]

# temp <- timegroups_names$timegroups[data_filt$timegroups]

cols_temp <- c("gridg1", "gridg3", "month", "timegroups")

data_2 <- data_filt %>% 
  mutate(across(.cols = all_of(cols_temp), ~ as.factor(.)))

rm(cols_temp)

# start parallel
# n.cores = worker_procs # From command line

# if(have_run_stats) {
#   run_stats <- species_run_stats
# }
species_todo <- length(species_to_process)
# if (species_todo==0) {
#   species_todo <- length(listofspecies)
# } else {
#   if(have_run_stats) {
#     run_stats <- run_stats %>%
#       filter(species_name %in% species_to_process)
#   }
# }
#message(paste("Processing", species_todo, "species..."))

# species_threads <- min(n.cores, species_todo) # if very few species then can't engage all cores
# species_threads_active <- 0
# species_done <- 0
# species_failed <- 0
trends0 <- NULL

# run_stats is ordered in descending order of runtime.
# longest running species typically consume max memory as well
# so peak memory consumption should happen in the beginning
# then this will taper.  Doing this also ensures that the job
# will fail in the beginning rather than later due to OOM
# scenarios.
#
# the table has species name, time to run, and peakRAM usage
# for that run.  Having all this as data makes it possible to
# "schedule" intelligently later

#free_ram <- get_free_ram() - ram_safety_margin*1000*1024
# 
# if(have_run_stats) {
#   message("rs = ", have_run_stats)
#   # Minimum we can run is 1 species at a time. If we don't have memory
#   # for that, better to exit now!
#   min_ram_needed <- max(run_stats$peakRAM)*1000*1024
#   
#   if (min_ram_needed > free_ram) {
#     message("Not enough RAM to fit even single species.")
#     message(paste("Min reqd:", min_ram_needed, "Available:", free_ram))
#     message("Increase container memory limits and try again.")
#     quit()
#   }
# }

#message("Free RAM at start ", as.integer(free_ram/1000000), " MB")
# if(!have_run_stats) {
#   message("No run stats. Running jobs in FCFS order.")
# } else if(ram_interleave) {
#   message("Running jobs with runtime length, combined with RAM interleave scheduling")
#   # bucket based on RAM. Mix of jobs with different RAM (3+2+1)=6 = 2x3 cores
#   # or 2 GB ram per core
#   bucket1 <- subset(run_stats, peakRAM < 1000)
#   bucket2 <- subset(run_stats, peakRAM >= 1000 & peakRAM < 2000)
#   bucket3 <- subset(run_stats, peakRAM >= 2000)
#   
#   # clear
#   run_stats <- run_stats[0, ]
#   # combine buckets with interleaving
#   iters <- max(nrow(bucket1),nrow(bucket2),nrow(bucket3))
#   for (i in 1:iters) {
#     if(i<=nrow(bucket3)) {
#       run_stats[nrow(run_stats)+1,] <- bucket3[i,]
#     }
#     if(i<=nrow(bucket2)) {
#       run_stats[nrow(run_stats)+1,] <- bucket2[i,]
#     }
#     if(i<=nrow(bucket1)) {
#       run_stats[nrow(run_stats)+1,] <- bucket1[i,]
#     }
#   }
# } else {
#   message("Running jobs with runtime length scheduling")
# }
# 
# try_thread_start <- TRUE
# 
# launched <- list()
# if(have_run_stats) {
#   species_pending_list <- as.vector(run_stats$species_name)
# } else {
#   if(length(species_to_process)==0) {
#     species_pending_list <- listofspecies
#   } else {
#     species_pending_list <- species_to_process
#   }
# }

# Keep going as long as we have species to run, or threads active
# start as many threads as we have (remaining) capacity for
started <- 0

# if(have_run_stats) {
#   launch_species_idx <- which(run_stats$species_name == launch_species)
#   min_ram_needed <- run_stats$peakRAM[launch_species_idx]*1000*1024
#   if(min_ram_needed > free_ram) {
#     message("Won't start ", species_threads-species_threads_active,
#             " threads due to insufficient RAM (needed for 1 more: ",
#             as.integer(min_ram_needed/1000000), " MB free:",
#             as.integer(free_ram/1000000), " MB)")
#     break
#   }
# } else {
#   min_ram_needed <- 0
# }
launch_species <- species_to_process[1]
species_index <- which(species_names$COMMON.NAME==species_to_process[1])

# assume job will consume this much
# free_ram <- free_ram - min_ram_needed

# proc <- singlespeciesrun(
#     container = container,
#     reproducible = reproducible_run,
#     stats_dir = trends_stats_dir,
#     species_dir = trends_species_dir,
#     data = data_2,
#     species_index = species_index,
#     species = launch_species,
#     specieslist = specieslist, 
#     restrictedspecieslist = restrictedspecieslist,
#     singleyear = singleyear)

container = container
reproducible = reproducible_run
#stats_dir = trends_stats_dir
#species_dir = trends_species_dir
data = data_2
species_index = species_index
species = launch_species
specieslist = specieslist
restrictedspecieslist = restrictedspecieslist
singleyear = singleyear

# retval <- singlespeciesrun_internal(container, reproducible, data, species_index, species,
#                                     specieslist, restrictedspecieslist, singleyear)


data1 = data # Data1 is filtered for the random group IDs.
rm(data)

# Checkpoint
data1_latest <- data1
save(data1_latest, file = "01_analyses_full/data1_before_filter_sp_latest.RData")

if(reproducible) {
  message("Setting seed to 0 to ensure reproducible runs")
  set.seed(0)
}

# get information for the species of interest 
specieslist2 = specieslist %>% filter(COMMON.NAME == species)

# three different flags for three different model types that will be run.
# 0 is normal model, with full random effects. depending on restricted species,
# model changes slightly.
flag = 0
# if (species %in% restrictedspecieslist$COMMON.NAME)
# {
#   flag = 1
#   restrictedlist1 = restrictedspecieslist %>% filter(COMMON.NAME == species)
#   specieslist2$ht = restrictedlist1$ht
#   specieslist2$rt = restrictedlist1$rt
#   
#   if (restrictedlist1$mixed == 0) {
#     flag = 2
#   }
# }

# filters data based on whether the species has been selected for long-term trends (ht) 
# or short-term trends (rt) 
# (if only recent, then need to filter for recent years. else, use all years so no filter.)

# if (singleyear == FALSE) {
#   
#   if (is.na(specieslist2$ht) & !is.na(specieslist2$rt)) {
#     data1 = data1 %>% filter(year >= soib_year_info("cat_start", container))
#   }
#   
# } else if (singleyear == TRUE) {
#   
#   data1 = data1 %>% filter(year == soib_year_info("latest_year", container))
# }


data1 = data1 %>%
  filter(COMMON.NAME == species_index) %>%
  distinct(gridg3, month) %>% 
  left_join(data1) %>%
  suppressMessages()

# Checkpoint
data1_grfr_list_info_latest <- data1
save(data1_grfr_list_info_latest, file = "01_analyses_full/data1_after_filter_sp_list_info_latest.RData")

data1_grfr_latest = data1 %>%
  filter(COMMON.NAME == species_index)

save(data1_grfr_latest, file = "01_analyses_full/data1_after_filter_sp_latest.RData")

dataset_size = nrow(data1)

colnames(data1)

tm = data1 %>% distinct(timegroups)
#rm(data, pos = ".GlobalEnv")

datay = data1 %>%
  distinct(gridg3, gridg1, group.id, .keep_all = TRUE) %>% 
  group_by(gridg3, gridg1) %>% 
  reframe(medianlla = median(no.sp)) %>%
  group_by(gridg3) %>% 
  reframe(medianlla = mean(medianlla)) %>%
  reframe(medianlla = round(mean(medianlla)))

medianlla = datay$medianlla

# expand dataframe to include absences as well

setDT(data1)

# Get distinct rows and filter based on a condition
# (using base data.table because lazy_dt with immutable == FALSE would
# modify data even though we are assigning to checklistinfo.
# and immutable == TRUE copies the data and this is a huge bottleneck)
# considers only complete lists

# Added OBSERVER.ID
checklistinfo <- unique(data1[, 
                              .(gridg1, gridg3, ALL.SPECIES.REPORTED, OBSERVER.ID,
                                group.id, month, year, no.sp, timegroups)
])[
  # filter
  ALL.SPECIES.REPORTED == 1
]

checklistinfo <- checklistinfo[
  , 
  .SD[1], # subset of data
  by = group.id
]

checklistinfo_latest <- checklistinfo

# Checkpoint

save(checklistinfo_latest, file = "01_analyses_full/checklistinfo_latest.RData")


# expand data frame to include the bird species in every list

join_by_temp <- c("group.id", "gridg1", "gridg3",
                  "ALL.SPECIES.REPORTED", "OBSERVER.ID", "month", "year",
                  "no.sp", "timegroups", "COMMON.NAME")

expanded_broken = checklistinfo %>% 
  lazy_dt(immutable = FALSE) |> 
  mutate(COMMON.NAME = species_index) %>% 
  left_join(data1 |> lazy_dt(immutable = FALSE),
            by = join_by_temp) %>%
  dplyr::select(-c("COMMON.NAME",
                   "ALL.SPECIES.REPORTED","group.id","year")) %>%
  # deal with NAs (column is character)
  mutate(OBSERVATION.COUNT = case_when(is.na(OBSERVATION.COUNT) ~ 0,
                                       OBSERVATION.COUNT != 0 ~ 1,
                                       TRUE ~ as.numeric(OBSERVATION.COUNT))) |> 
  as_tibble()

expanded_latest_broken <- expanded_broken

# Checkpoint

save(expanded_latest_broken, file = "01_analyses_full/expanded_latest_broken.RData")

rm(join_by_temp)


ed_broken = expanded_broken %>%
  # converting months to seasons
  mutate(month = as.numeric(month)) %>% 
  mutate(month = case_when(month %in% c(12,1,2) ~ "Win",
                           month %in% c(3,4,5) ~ "Sum",
                           month %in% c(6,7,8) ~ "Mon",
                           month %in% c(9,10,11) ~ "Aut")) %>% 
  mutate(month = as.factor(month))

colnames(ed_broken)

# save some values referenced later so we can get rid of memory hog data1
gg1 <- data1$gridg1[1]
gg3 <- data1$gridg3[1]
rm(data1)

# the model ---------------------------------------------------------------

fixed_effects <- "OBSERVATION.COUNT ~ month + month:log(no.sp)"
include_timegroups <- if (singleyear == FALSE) "+ timegroups" else 
  if (singleyear == TRUE) ""
random_effects <- if (flag == 0) "+ (1|gridg3/gridg1)" else 
  if (flag == 1) "+ (1|gridg1)" else 
    if (flag == 2) ""

model_formula <- as.formula(glue("{fixed_effects} {include_timegroups} {random_effects}"))

m1_broken <- 
  glmer(model_formula, 
        data = ed_broken, family = binomial(link = 'cloglog'), 
        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))



# predicting from model ---------------------------------------------------

# prepare a new data file to predict

ltemp_broken <- ed_broken %>% 
  {if (singleyear == FALSE) {
    group_by(., month) %>% 
      reframe(., timegroups = unique(tm$timegroups))
  } else if (singleyear == TRUE) {
    distinct(., month)
  }} %>% 
  mutate(no.sp = medianlla,
         # taking the first value but any random value will do because we do not
         # intend to predict random variation across grids
         gridg1 = gg1,
         gridg3 = gg3)

f2_broken <- ltemp_broken %>% 
  {if (singleyear == FALSE) {
    dplyr::select(., timegroups)
  } else if (singleyear == TRUE) {
    .
  }} %>% 
  mutate(freq = 0, se = 0) %>%  # this is not actually needed
  {if (singleyear == FALSE) {
    .
  } else if (singleyear == TRUE) {
    dplyr::select(., freq, se)
  }}

  #pred = predict(m1, newdata = ltemp, type = "response", re.form = NA, allow.new.levels=TRUE)
  pred_broken = predictInterval(m1_broken, newdata = ltemp_broken, which = "fixed",
                         level = 0.48, type = "linear.prediction",
                         include.resid.var = TRUE)
  f2_broken$freqt = pred_broken$fit
  f2_broken$set = pred_broken$fit-pred_broken$lwr


f1_broken = f2_broken %>%
  filter(!is.na(freqt) & !is.na(set)) %>%
  # average across month
  {if (singleyear == FALSE) {
    group_by(., timegroups) %>% 
      reframe(freq = mean(freqt), se = mean(set)) %>% 
      right_join(tm) %>% 
      left_join(databins %>% distinct(timegroups, year)) %>% 
      rename(timegroupsf = timegroups,
             timegroups = year) %>% 
      mutate(timegroupsf = factor(timegroupsf, 
                                  levels = soib_year_info("timegroup_lab", container))) %>%
      complete(timegroupsf) %>% 
      arrange(timegroupsf) %>%
      suppressMessages()
  } else if (singleyear == TRUE) {
    reframe(., freq = mean(freqt), se = mean(set))
  }}



tocomb_broken = c(species, f1_broken$freq, f1_broken$se)

# keep track of which PID is which species
# if(have_run_stats) {
#   run_duration <-run_stats$time[launch_species_idx]
# } else {
#   run_duration <- 0
# }
# launched[[as.character(proc$pid)]] <- c(launch_species,
#                                         proc.time()[3],
#                                         run_duration)
# 
# if(have_run_stats) {
#   message("Started: ",
#           launch_species, ", estimated peakRAM: ",
#           as.integer(min_ram_needed/1000000), " MB,",
#           " runtime: ", run_duration, " seconds")
# } else {
#   message("Started: ", launch_species)
# }

# species_threads_active <- species_threads_active + 1
# started <- started + 1


#      if(started > 0) {
#        message(paste("Threads started:", started))
#      }

# if(have_run_stats) {
#   message("T=",proc.time()[3],
#           " Estmated Peak Free RAM:", as.integer(free_ram/1000000), " MB",
#           " Threads:", species_threads_active,
#           " Done:", species_done,
#           " Pending:", length(species_pending_list),
#           " Failed:", species_failed)
# } else {
#   message("T=",proc.time()[3],
#           " Threads:", species_threads_active,
#           " Done:", species_done,
#           " Pending:", length(species_pending_list),
#           " Failed:", species_failed)
# }
# 
# done <- 0
# result_set <- mccollect(timeout = 1000, wait=FALSE)
# if (!is.null(result_set)) {
#   finished_pids <- names(result_set)
#   # iterate over the results
#   for (idx in 1:length(finished_pids)) {
#     result <- result_set[[finished_pids[idx]]]
#     pid_str <- as.character(finished_pids[idx])
#     this_species <- launched[[pid_str]][1]
#     time_taken <- proc.time()[3] - as.numeric(launched[[pid_str]][2])
#     time_expected <- as.numeric(launched[[pid_str]][3])
#     launched[[pid_str]] <- NULL # remove
#     if(length(result)==0) {
#       message("Failed: ", this_species, " Time taken:", time_taken, " secs")
#       # Move this failed to pending
#       species_pending_list <- append(species_pending_list, this_species)
#       species_failed <- species_failed + 1
#     } else {
#       # show 100% without run stats
#       if(!have_run_stats) {
#         time_expected <- time_taken
#       }
#       percent <- (time_taken/time_expected)*100.0
#       percent_str <- format(round(percent, 2), nsmall = 2)
#       message("Finished: ", this_species, " Time taken:", time_taken, " secs (", percent_str, " %)")
#       trends0 <- cbind(trends0, result)
#       species_done <- species_done + 1
#     }
#     
#     if(have_run_stats) {
#       done_index <- which(run_stats$species_name==this_species)
#       done_mem <- run_stats$peakRAM[done_index]*1000*1024
#       free_ram <- free_ram + done_mem
#     }
#     species_threads_active <- species_threads_active - 1
#     done <- done + 1
#   }
#   # Thread status changed, so may need to launch next
#   try_thread_start <- TRUE
# } else {
#   try_thread_start <- FALSE
# }
# 
# #      if (done > 0) {
# #        message(paste("Threads finished:", done))
# #      }



trends = data.frame(tocomb_broken) %>% 
  # converting first row of species names (always true) to column names
  magrittr::set_colnames(.[1,]) %>% 
  slice(-1) %>% 
    mutate(.,
           timegroupsf = rep(databins$timegroups, 2),
           timegroups = rep(databins$year, 2),
           type = rep(c("freq", "se"), 
                      # will always have 2*N.YEAR rows (freq, se)
                      each = length(soib_year_info("timegroup_lab", container))),
           sl = k) %>%  # sim number
      # pivoting species names longer
      pivot_longer(-c(timegroups, timegroupsf, sl, type), 
                   values_to = "value", names_to = "COMMON.NAME") %>% 
      pivot_wider(names_from = type, values_from = value) %>% 
      # numerical ID for species names, for arranging
      mutate(sp = row_number(), .by = timegroupsf) %>%
      arrange(sl, sp) %>%
      dplyr::select(-sp) %>% 
      # reordering
      relocate(sl, COMMON.NAME, timegroupsf, timegroups, freq, se) %>% 
  # make sure freq and se are numerical
  mutate(across(c("freq", "se"), ~ as.numeric(.)))

# reorder species to reflect the list order prior to RAM optimizations
# this makes it easy to compare trends_N.csv files with older ones side
# by side
# if full run, overwrite the CSV
# else append single year results to all previous year results
  
  write.csv(trends, file = write_path, row.names = FALSE)
  



tictoc::toc() 

# maybe unnecessary time-consuming step:
# worth trying the profiling mentioned here http://adv-r.had.co.nz/memory.html
# https://stackoverflow.com/questions/1467201/forcing-garbage-collection-to-run-in-r-with-the-gc-command
# gc()


