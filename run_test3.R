day_number = 2
day_number_base = as.character(day_number)
day_base = "DAY"
day_dir = paste(day_base,day_number_base,sep='')
if (!dir.exists(day_dir))
  dir.create(day_dir)
# setwd(day_dir)

