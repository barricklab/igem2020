
# This function calculates the burst size of a given simulation 
burst_size_function <- function(sim_data, file_name, lysis_time){
  
  # read in file that defines the number of proteins (copies) required per gene per virion
  vreq = read_csv("proteins_in_t7_virion.csv")
  
  #pick genes for virion proteins (only works for output from latest script update since gene labels were changed to gp)
  virion.species = c("gp6.7", "gp7.3", "gp8", "gp10A", "gp11", "gp12", "gp14", "gp15", "gp16", "gp17")
  
  # checks if lysis_time is -1 (does not lyse), if it is, set lysis time to the end of simulation
  if (lysis_time == -1){
    lysis_time_num = max(sim_data$time)
  }
  else{
    lysis_time_num = lysis_time
  }
  
  # Creates datatset filtering for the lysis time and the species listed above
  virion_genes = sim_data %>%  
    select(species, protein, time) %>%
    filter(species %in% virion.species) %>%
    filter(time <= lysis_time_num)
  
  # selects max protein count for each virion gene 
  vt <- sqldf("select max(protein), species from virion_genes group by species")
  
  # joins the two input files by the mutual species type column, in order to make the following calculations. 
  virion.table <- left_join(vt, vreq,by="species")
  
  # creates a new column that calculates the max number of virions produced per protein type (in reference to the t7 virion file)
  virion.table$max_virions <- as.numeric(as.character(virion.table$`max(protein)`)) / virion.table$copies
  
  #finds the smallest number of total virions generated, which would satisfy all protein count requirements   
  total_virions_generated = min(virion.table[,"max_virions"])
  
  if (nrow(virion.table) != 10){
    total_virions_generated = 0
  }
  
  # gets filename prefix
  file_prefix = gsub(pattern = "\\.counts.tsv$", "", file_name)
  
  # generates a dataframe storing the file prefix and the virions generated
  test.data <- data.frame(file_prefix, total_virions_generated, lysis_time)
  
  
  return(test.data)
  
}


########Function to calculate lysis time

find_lts <- function(sim_data){
  
  #create df with only holin species, create df with only lysozyme species
  mt_holin <- filter(sim_data, species %in% c("gp17.5"))
  mt_lysin <- filter(sim_data, species %in% c("gp3.5"))
  
  
  #extract rows where holin protein > protein required for lysis, same with lysozyme
  holin_passed <- mt_holin %>% filter(protein>=12341)
  
  
  #if there is data in holin_passed & lysin_passed, find the smallest time -> this is a potential lysis time
  if (nrow(holin_passed)>=1 && nrow(mt_lysin)>=1 ){
    lt_holin <- min(holin_passed$time)
    lt_final <- lt_holin
  }
  
  #if there is no data in holin_passed, then holin levels either
  #1. have not reached required level for lysis or
  #2. holin gene was deleted
  #this is testing for case #1: holin is present but hasn't reached critical conc. yet -> longer sim. required.
  if (nrow(holin_passed)==0 && nrow(mt_holin)>=1){
    lt_final <- -1
    
  }
  
  #if holin gene was deleted, lysis can still occur, potentially through gp19.5
  
  
  #if lysin gene was deleted, lysis can still occur, potentially through gp16
  #NOT BEEN UPDATED BY TRANSCRIPT DEGRADATION
  if (nrow(mt_lysin)==0){
    gp16 <- sim_data %>% filter (species %in% c("gp16"), protein>=19744)
    if (nrow(gp16)>=1){
      lt_final <- min(gp16$time)
    }
    else {
      lt_final <- -1
    }
  }
  
  #if gene 19.5 is deleted, lysis is delayed
  
  #return final lysis time -> will be either time or -1
  return(lt_final)
  
}

