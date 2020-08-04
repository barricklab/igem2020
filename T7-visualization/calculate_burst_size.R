require(sqldf)
require(dplyr)
require(readr)

# This function calculates the burst size of a given simulation 
burst_size_function <- function(file_name, lysis_time){
  
  # reads in count.tsv file900
  all_data = read_tsv(file_name)
  
  # read in file that defines the number of proteins (copies) required per gene per virion
  vreq = read_csv("proteins_in_t7_virion.csv")
  
  #pick genes for virion proteins (only works for output from latest script update since gene labels were changed to gp)
  virion.species = c("gp6.7", "gp7.3", "gp8", "gp10A", "gp11", "gp12", "gp14", "gp15", "gp16", "gp17")
  
  # Creates datatset filtering for the lysis time and the species listed above
  virion_genes = all_data %>%  
    select(species, protein, time) %>%
    filter(species %in% virion.species) %>%
    filter(time <= lysis_time)
  
  # selects max protein count for each virion gene 
  vt <- sqldf("select max(protein), species from virion_genes group by species")
  
  # joins the two input files by the mutual species type column, in order to make the following calculations. 
  virion.table <- left_join(vt, vreq,by="species")
  
  # creates a new column that calculates the max number of virions produced per protein type (in reference to the t7 virion file)
  virion.table$max_virions <- as.numeric(as.character(virion.table$`max(protein)`)) / virion.table$copies
  
  #finds the smallest number of total virions generated, which would satisfy all protein count requirements   
  total_virions_generated = min(virion.table[,"max_virions"])
  
  # gets filename prefix
  file_prefix = gsub(pattern = "\\.counts.tsv$", "", file_name)
  
  # generates a dataframe storing the file prefix and the virions generated
  test.data <- data.frame(file_prefix, total_virions_generated, lysis_time)

  
  return(test.data)
  
}

# Creates a list with all the .counts.tsv files in the current working directory. If you want to get files from
# another directory other than current directory, type out the path inside the parenthesis of list.files()
all.files <- list.files()
tsv.files <- grep(".counts.tsv",all.files,value=T)

# Initialize dataframe to store burst size from all the files in directory
virion.data = data.frame()

# asks user for the lysis time (applies to all files inputed). For now, just enter 1500.
lysis_time <- as.numeric(readline("Enter the lysis time: "))


# loops through each file and calls function defined above and appends test.data to virion.data
for(i in tsv.files){
  # eventually may want to call lysis time function here to set lysis time for each individual file?
  test.data = burst_size_function(i, lysis_time)
  virion.data = rbind(virion.data, test.data)
  
}

# creates csv with virion.data and saves it in current directory. 
# change the file name depending on what data you are running. 
write_csv(virion.data, "WT_Burst_Size_1500s.csv")







