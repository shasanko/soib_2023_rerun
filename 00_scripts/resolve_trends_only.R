library(tictoc)
library(VGAM)
library(sf)
library(writexl)
library(tidyverse)

source('00_scripts/00_functions.R')
cur_mask <- "none"
interannual_update <- FALSE
run_res_trends <- TRUE

species_to_process <- c("Gray Francolin",
                        "Ashy Prinia",
                        "Black-winged Kite")

# setup -------------------------------------------------------------------

# preparing data for specific mask (this is the only part that changes, but automatically)
cur_metadata <- get_metadata(cur_mask)

# read paths
base_path <- cur_metadata$FULLSPECLIST.PATH
speclist_path <- cur_metadata$SPECLISTDATA.PATH
trends_pathonly <- cur_metadata$TRENDS.PATHONLY

# write paths
lttsens_path <- cur_metadata$LTTSENS.PATH
cursens_path <- cur_metadata$CURSENS.PATH
trends_outpath <- cur_metadata$TRENDS.OUTPATH

mainwocats_path <- cur_metadata$SOIBMAIN.WOCATS.PATH

recentcutoff = soib_year_info("cat_start")

load("00_data/spec_misid.RData") # to remove from LTT and CAT "selection" later
# for occupancy
load("00_data/maps_sf.RData")
load(speclist_path)

# data processing and prep ------------------------------------------------

base = read.csv(base_path) %>% 
  filter(COMMON.NAME %in% species_to_process) %>%
  # if full column has no X at all, gets read as NAs
  mutate(across(c(Long.Term.Analysis, Current.Analysis, Selected.SOIB),
                ~ as.character(.))) %>%
  mutate(across(c(Long.Term.Analysis, Current.Analysis, Selected.SOIB),
                ~ replace_na(., ""))) %>%
  dplyr::select(-SCIENTIFIC.NAME)

base = base %>% dplyr::select(-c("X"))

main = read.csv("00_data/SoIB_mapping_2022.csv") %>%
  filter(eBird.English.Name.2022 %in% species_to_process) %>%
  left_join(base, by = c("eBird.English.Name.2022" = "COMMON.NAME"))

# separate object for sensitivity analysis
sens = main %>% dplyr::select(eBird.English.Name.2022)

# trends files

trends <- list.files(path = "~/new_pipeline/soib_2023/01_analyses_full/trends/", 
                     # Generate the full file paths
                     full.names = T) 

# The output contains both folders (species_1, species_2....) and files. We
# just need the files

trends_files <- trends[!file.info(trends)$isdir]


trends <- trends_files %>% 
  # Read each CSV file and combine them into a single data frame
  map_df(read.csv)

# data filtering: problem species -----------------------------------------

# rewriting selected species for LTT and CAT
spec_lt = main %>%
  filter(Long.Term.Analysis == "X") %>%
  pull(eBird.English.Name.2022)

spec_ct = main %>%
  filter(Current.Analysis == "X") %>%
  pull(eBird.English.Name.2022)


# checkpoint-object "main"
main1_postfilt <- main


# calculations: prep --------------------------------------------------------

# This section performs calculations on the trends data frame to derive new columns, 
# such as lci, mean, and rci. It also creates a data frame called trends_framework to 
# define the timegroups and species combinations. (Summary from ChatGPT)

# Mean is calculated first, then CI. This is done for long-term and current 
# trend separately.


# Years to project for PJ's IUCN comparison
extra.years = soib_year_info("iucn_projection")

# trends_reju <- trends %>% 
#   filter(COMMON.NAME == "Red Junglefowl" & timegroupsf == "2000-2006")
# mean(trends_reju$freq)

trends_combined = trends %>%
  group_by(COMMON.NAME, timegroupsf,  timegroups) %>%
  reframe(mean_trans = mean(freq))

trends = trends %>%
  group_by(COMMON.NAME, timegroupsf, timegroups) %>% 
  reframe(mean_trans = mean(freq), 
          # adding se from variation between the means and propagated SE
          se_trans = sd(freq) + sqrt(sum(se^2)/n())) %>%
  group_by(COMMON.NAME, timegroupsf, timegroups) %>% 
  mutate(lci = clogloglink(mean_trans - 1.96*se_trans, inverse = T),
         mean = clogloglink(mean_trans, inverse = T),
         rci = clogloglink(mean_trans + 1.96*se_trans, inverse = T)) %>%
  ungroup()

# adding extra years to dataframe (and completing species combinations)
trends_framework <- trends %>%
  distinct(timegroups, COMMON.NAME) %>%
  tidyr::expand(timegroups = c(unique(timegroups), extra.years),
                COMMON.NAME) %>%
  complete(timegroups, COMMON.NAME) %>% 
  left_join(trends %>% distinct(COMMON.NAME, timegroupsf, timegroups)) %>% 
  mutate(timegroupsf = ifelse(is.na(timegroupsf), timegroups, timegroupsf))


modtrends = na.omit(trends) %>% # NAs are all spp. not included in long-term
  filter(COMMON.NAME %in% spec_lt) %>%
  arrange(COMMON.NAME, timegroups) %>%
  # getting trends values of first year (only for long-term, i.e., pre-2000)
  group_by(COMMON.NAME) %>% 
  # _trans are link-scale, "mean" is back-transformed
  mutate(m1 = first(mean_trans),
         mean_year1 = first(mean),
         s1 = first(se_trans)) %>% 
  ungroup() %>% 
  # for calculating change in abundance index (as % change)
  mutate(mean_std = 100*mean/mean_year1) # back-transformed so value is % of year1 value

# sensitivity check to ensure edge species are later converted to the conservative status
modtrends1 <- ltt_sens_sim(my_seed = 1, data = modtrends)
modtrends2 <- ltt_sens_sim(my_seed = 2, data = modtrends)
modtrends3 <- ltt_sens_sim(my_seed = 3, data = modtrends)
modtrends4 <- ltt_sens_sim(my_seed = 4, data = modtrends)
modtrends5 <- ltt_sens_sim(my_seed = 5, data = modtrends)
save(modtrends1, modtrends2, modtrends3, modtrends4, modtrends5, file = lttsens_path)

# "main" simulated CIs
set.seed(10) 
modtrends = modtrends %>% 
  # calculating CIs
  group_by(COMMON.NAME, timegroups) %>% 
  # 1000 simulations of transformed ratio of present:original values
  # quantiles*100 from these gives us our CI limits for mean_std
  reframe(tp0 = simerrordiv(mean_trans, m1, se_trans, s1)$rat) %>% 
  group_by(COMMON.NAME, timegroups) %>% 
  reframe(lci_std = 100*as.numeric(quantile(tp0, 0.025)),
          rci_std = 100*as.numeric(quantile(tp0, 0.975))) %>% 
  right_join(modtrends, by = c("COMMON.NAME", "timegroups"))

# saving the values for final year in "main" as well:
# temp object then left_join instead of right_join because species order in main
# needs to be preserved
temp <- modtrends %>%
  filter(timegroups == soib_year_info("latest_year")) %>%
  dplyr::select(COMMON.NAME, lci_std, mean_std, rci_std) %>%
  rename(longtermlci = lci_std,
         longtermmean = mean_std,
         longtermrci = rci_std)

main <- main %>%
  left_join(temp, by = c("eBird.English.Name.2022" = "COMMON.NAME"))

modtrends = modtrends %>%
  group_by(COMMON.NAME) %>%
  # making CI band zero for first year
  mutate(lci_std = case_when(timegroups == first(timegroups) ~ mean_std,
                             TRUE ~ lci_std),
         rci_std = case_when(timegroups == first(timegroups) ~ mean_std,
                             TRUE ~ rci_std)) %>%
  ungroup() %>%
  dplyr::select(timegroupsf, timegroups, COMMON.NAME, lci_std, mean_std, rci_std)


# checkpoint-object "main"
main2_postLTT <- main

#}

# calculations: trends (current) ------------------------------------------

# Unlike trends for long-term, this is not last year divided by first year, instead a slope.
# For each year, sim 1000 points within CI for each year. In each sim, line is fitted.
# In each case, three models fitted (main, exponential for projection, sensitivity).
# Exponential is for Praveen's IUCN estimation using annual change, not used here.

tic("Calculating current trends")

# First, getting mean + CI for each year
modtrends_recent = trends %>%
  filter(COMMON.NAME %in% spec_ct &
           timegroups >= recentcutoff) %>%
  arrange(COMMON.NAME, timegroups) %>%
  # getting trends values of first year (only for current)
  group_by(COMMON.NAME) %>%
  mutate(m1 = first(mean_trans),
         mean_year1 = first(mean),
         s1 = first(se_trans)) %>%
  ungroup() %>%
  # for calculating change in abundance index (as % change)
  mutate(mean_std_recent = 100*mean/mean_year1)

# Getting CIs for each year
set.seed(10)
modtrends_recent = modtrends_recent %>%
  # calculating CIs
  group_by(COMMON.NAME, timegroups) %>%
  # 1000 simulations of transformed ratio of present:original values
  # quantiles*100 from these gives us our CI limits for mean_std
  reframe(simerrordiv(mean_trans, m1, se_trans, s1)) %>%  # gives rat and val columns
  group_by(COMMON.NAME, timegroups) %>%
  reframe(lci_std_recent = 100*as.numeric(quantile(rat, 0.025)),
          rci_std_recent = 100*as.numeric(quantile(rat, 0.975)),
          rat = rat,
          val = val) %>%
  right_join(modtrends_recent, by = c("COMMON.NAME", "timegroups"))

# Now running simulations for:
#   - slope and SE of main trend
#   - 8x slopes and SEs for sensitivity analyses
#   - exponential model predictions into future years

set.seed(1)
temp_sims <- modtrends_recent %>%
  group_by(COMMON.NAME) %>%
  dplyr::select(timegroups, rat, val) %>%
  group_by(COMMON.NAME, timegroups) %>%
  mutate(sim = 1:n(),
         # for each sim, sampling from 1000 val values within species-year group
         val_sample = map_dbl(sim, ~ sample(val, 1)))


# IN NEXT SECTION:
# Take these simulated values of reporting frequency for each year,
# fit lm through them, then predict new values based on that lm.

# Numerator and denominator for slope calculation, along with error propagation
# Uncertainty between 2015 and 2016 is highest, hence indexing them, i.e.,
# using that uncertainty for entire slope (being cautious).

# errordiv calculates and returns slope and SE from predicted values.


# making CI band zero for first year
modtrends_recent = modtrends_recent %>%
  # now we don't need rat and val
  distinct(COMMON.NAME, timegroupsf, timegroups,
           lci_std_recent, mean_std_recent, rci_std_recent) %>%
  group_by(COMMON.NAME) %>%
  mutate(lci_std_recent = case_when(timegroups == first(timegroups) ~ mean_std_recent,
                                    TRUE ~ lci_std_recent)) %>%
  mutate(rci_std_recent = case_when(timegroups == first(timegroups) ~ mean_std_recent,
                                    TRUE ~ rci_std_recent)) %>%
  ungroup() %>%
  dplyr::select(timegroupsf, timegroups, COMMON.NAME,
                lci_std_recent, mean_std_recent, rci_std_recent)


# calculations: trends (current): MAIN ------------------------------------

# At this point, we jave modtrends, modtrends_recent, temp_simps
# For testing, subset these 

sl_data_main <- temp_sims %>%
  group_by(COMMON.NAME, sim) %>%
  # group_modify() performs function per grouping, so sl and slse calculated for all sim
  group_modify(~ {
    
    datatopred <- .x %>% dplyr::select(timegroups)
    
    modelfit <- lm(val_sample ~ timegroups, data = .x)
    pred <- predict(modelfit, newdata = datatopred, se = TRUE)
    
    num <- pred$fit[2] - pred$fit[1]
    den <- abs(pred$fit[1])
    numse <- sqrt(pred$se.fit[1]^2 + pred$se.fit[2]^2)
    dense <- pred$se.fit[1]
    
    # indexing for mean and se
    sl <- 100 * errordiv(num, den, numse, dense)[1] %>% as.numeric()
    slse <- errordiv(num, den, numse, dense)[2] %>% as.numeric()
    
    .x %>%
      reframe(sl = sl,
              slse = slse)
    
  }) %>%
  group_by(COMMON.NAME) %>%
  reframe(mean.slope = mean(sl),
          se.slope = sd(sl) + sqrt(sum(slse^2)/length(slse))) %>%
  group_by(COMMON.NAME) %>%
  reframe(currentslopelci = mean.slope - 1.96*se.slope,
          currentslopemean = mean.slope,
          currentsloperci = mean.slope + 1.96*se.slope)

# joining to main data
main <- main %>%
  left_join(sl_data_main, by = c("eBird.English.Name.2022" = "COMMON.NAME"))


# checkpoint-object "main"
main3_postCATmain <- main

# calculations: trends (current): SENS ------------------------------------

# here, fitting the same linear model as for main trend, but for data in which
# one year is dropped each time.
# this is to see how much each year affects the estimate.
sl_data_sens <- temp_sims %>%
  group_by(COMMON.NAME, sim) %>%
  group_modify(~ {
    
    temp_df <- .x # to refer inside the second purrr function below
    
    .x %>%
      reframe(
        
        # produce as many iterations as no. of CAT years (that many sl and slse cols)
        !!!imap(soib_year_info("cat_years"), ~ {
          
          modelfit <- lm(val_sample ~ timegroups,
                         # one year dropped
                         data = temp_df[temp_df$timegroups != soib_year_info("cat_years")[.y],])
          
          pred <- predict(
            modelfit, se = TRUE,
            # one year dropped
            newdata = data.frame(
              timegroups = temp_df$timegroups[temp_df$timegroups != soib_year_info("cat_years")[.y]]
            )
          )
          
          num <- pred$fit[2] - pred$fit[1]
          den <- abs(pred$fit[1])
          numse <- sqrt(pred$se.fit[1]^2 + pred$se.fit[2]^2)
          dense <- pred$se.fit[1]
          
          # create col names based on iteration index
          col_name_sl <- paste0("sl", .y)
          col_name_slse <- paste0("slse", .y)
          
          tibble(
            !!col_name_sl := 100 * errordiv(num, den, numse, dense)[1] %>% as.numeric(),
            !!col_name_slse := errordiv(num, den, numse, dense)[2] %>% as.numeric()
          )
          
        }) %>% 
          bind_cols()
        
      )
    
  }) %>%
  
  # mean and SE of slope
  group_by(COMMON.NAME) %>%
  reframe(across(
    
    .cols = matches("^sl\\d+$"), # "sl"s but not "slse"s
    .fn = list(mean.slope = ~ mean(.),
               se.slope = ~ sd(.) + sqrt(sum(get(str_replace(cur_column(), "sl", "slse"))^2) /
                                           length(.))),
    .names = c("{.fn}_{.col}")
    
  )) %>%
  rename_with(~ str_remove(., "_sl")) %>%
  
  # (mean + SE) to (LCI, mean, RCI)
  group_by(COMMON.NAME) %>%
  reframe(across(
    
    .cols = starts_with("mean.slope"),
    .fn = list(lci = ~ . - 1.96 * get(str_replace(cur_column(), "mean", "se")),
               mean = ~ .,
               rci = ~ . + 1.96 * get(str_replace(cur_column(), "mean", "se"))),
    .names = c("currentslope{.fn}{.col}")
    
  )) %>%
  rename_with(~ str_remove(., "mean.slope"))

# joining to data object
sens <- sens %>%
  left_join(sl_data_sens, by = c("eBird.English.Name.2022" = "COMMON.NAME"))

write.csv(sens, 
          file = "01_analyses_full/results/current_sensitivity_1000.csv", 
          row.names = F)

# calculations: trends (current): PROJ ------------------------------------

ext_trends <- temp_sims %>%
  group_by(COMMON.NAME, sim) %>%
  group_modify(~ {
    
    datatopred <- data.frame(timegroups = extra.years)
    
    modelfit <- lm(log(val_sample) ~ timegroups, data = .x)
    pred <- predict(modelfit, newdata = datatopred, se = TRUE)
    
    # No need to calculate slope here.
    
    first_year <- .x %>% filter(timegroups == recentcutoff)
    
    .x %>%
      reframe(timegroups = datatopred$timegroups,
              mean = pred$fit,
              se = pred$se.fit,
              # value in first year
              val = first_year$val)
    
  }) %>%
  # we want mean and SE of change in repfreq between years, so divide by first year value
  mutate(lci_bt = 100*exp(mean - 1.96 * se)/val,
         mean_bt = 100*exp(mean)/val,
         rci_bt = 100*exp(mean + 1.96 * se)/val) %>%
  group_by(COMMON.NAME, timegroups) %>%
  reframe(lci_ext_std = mean(lci_bt),
          mean_ext_std = mean(mean_bt),
          rci_ext_std = mean(rci_bt))

toc()

#}

# calculations: combining two trends --------------------------------------

trends = trends %>%
  left_join(modtrends) %>%
  left_join(modtrends_recent) %>%
  dplyr::select(timegroups, COMMON.NAME, timegroupsf, mean_trans, se_trans,
                lci, mean, rci, lci_std, mean_std, rci_std,
                lci_std_recent, mean_std_recent, rci_std_recent) %>%
  right_join(trends_framework) %>%
  # projected trends
  left_join(ext_trends) %>%
  # changing names to factor, so they maintain order
  mutate(COMMON.NAME = factor(COMMON.NAME, levels = base$COMMON.NAME)) %>%
  arrange(COMMON.NAME, timegroups) %>%
  # replacing NAs
  mutate(lci_comb_std = case_when(is.na(lci_ext_std) ~ lci_std_recent,
                                  !is.na(lci_ext_std) ~ lci_ext_std),
         mean_comb_std = case_when(is.na(mean_ext_std) ~ mean_std_recent,
                                   !is.na(mean_ext_std) ~ mean_ext_std),
         rci_comb_std = case_when(is.na(rci_ext_std) ~ rci_std_recent,
                                  !is.na(rci_ext_std) ~ rci_ext_std)) %>%
  # truncating LCI at 0
  mutate(lci_comb_std = case_when(lci_comb_std < 0 ~ 0,
                                  TRUE ~ lci_comb_std)) %>%
  # ensuring correct order of columns
  relocate(timegroups, COMMON.NAME, timegroupsf, mean_trans, se_trans,
           lci, mean, rci, lci_std, mean_std, rci_std,
           lci_std_recent, mean_std_recent, rci_std_recent,
           lci_ext_std, mean_ext_std, rci_ext_std,
           lci_comb_std, mean_comb_std, rci_comb_std)


write.csv(trends, 
          file = "01_analyses_full/results/trends.csv", 
          row.names = F)

