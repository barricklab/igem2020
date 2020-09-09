require(sqldf)
require(dplyr)
require(readr)
require(tidyverse)
source("bs_lt_functions.R")

# Creates a list with all the .counts.tsv files in the current working directory. If you want to get files from
# another directory other than current directory, type out the path inside the parenthesis of list.files()
all.files <- list.files()
tsv.files <- grep(".counts.tsv",all.files,value=T)

# Initialize dataframe to store burst size from all the files in directory
model.data = data.frame()

# loops through each file, reads it, and calls function defined above and appends test.data to virion.data
for(i in tsv.files){
  file_name = i
  sim_data = read_tsv(i)
  lysis_time = find_lts(sim_data)
  test.data = burst_size_function(sim_data, file_name, lysis_time)
  model.data = rbind(model.data, test.data)
  
}

# creates csv with virion.data and saves it in current directory. 
# change the file name depending on what data you are running. 
write_csv(model.data, "WT_test.bslt.csv")


