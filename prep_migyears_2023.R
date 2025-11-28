
full_soib_my <- c(full_soib_my[1:51])
latest_soib_my <- 2022

save(full_soib_my, latest_soib_my, median_soib_hist_years, 
     file = "00_data/current_soib_migyears.RData")