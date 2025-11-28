# Remaps randomgroupids to 1000 separate files
# with integer mapped random group ids.
# Each run needs only one file
library("tictoc")

print("loading random group ids...")
tic("time")
load("01_analyses_full/randomgroupids.RData")
toc()

convert_group_id <- function(x) {
  base <- ifelse(substr(x,1,1)=="S",0L,1000000000L) # G=>giga base
  return(base + as.integer(substr(x,2,12)))
}

rgids <- randomgroupids
# for (i in 1:1000) {
#   remap <- convert_group_id(rgids[,i])
#   randomgroupids <- data.frame(group.id = remap)
#   print(paste("Random group iteration:", i, "of 1000"))
#   save(randomgroupids, file=paste0("01_analyses_full/rgids-",i,".RData"))
# }

for (i in 1:1000) {
  group.id <- convert_group_id(rgids[,i])
  randomgroupids <- group.id
  print(paste("Random group iteration:", i, "of 1000"))
  save(randomgroupids, file=paste0("01_analyses_full/rgids-",i,".RData"))
}

print("Done")
