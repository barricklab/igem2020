require (ggplot2)
require(readr)
source ("bs_lt_functions.R")

#creates list with all .counts.tsv files in working directory
runs = list.files(pattern = "*.tsv")

#initialize df for lysis time of each run, initialize a null dataframe for the nulls
lts = data.frame()
nulls = data.frame()

#runs the lyis time calculator through each file, adds lysis time to lts df
#if NULL (cell did not lyse during simulation time), add to nulls df
for (i in runs){
  file_name = i
  mt = read_tsv(file_name)
  lt = find_lts(mt)
  if (lt == -1)
  {
    nulls = rbind(nulls, "none")
  }
  else{
    lts = rbind(lts, lt)
  }
}
colnames(lts)<- c("lysis time")
mean_lt = mean(lts$`lysis time`)
median_lt = median(lts$`lysis time`)
sd_lt = sd(lts$`lysis time`)

nolyse = paste("Did not lyse: ", as.character(nrow(nulls)))
lyse = paste("Lysed: ", as.character(nrow(lts)))
mean_string = paste("Mean: ", as.character(mean_lt))
median_string = paste("Median: ", as.character(median_lt))
standarddeviation_string = paste("Std Dev.: ", as.character(sd_lt))


p <- hist(lts$`lysis time`,
          main = "Distribution of Lysis Time",
          xlab="Lysis Time")
p
title(sub=nolyse, adj = 1, line = 3, font = 2)
title(sub=lyse, adj=1, line = 2, font = 2)
mtext(mean_string, side = 3, line = 0, font = 0.5)
mtext(median_string, side =3, line = 0.7, font = 1)
mtext(standarddeviation_string, side =3, line = -.7, font = 1 )