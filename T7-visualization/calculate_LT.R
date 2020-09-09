require(dplyr)
require(tidyverse)

########Function to calculate lysis time

find_lts <- function(mutant){
  
  
  #read in mutant file
  mt = read_tsv(mutant)
  
  #create df with only holin species, create df with only lysozyme species
  mt_holin <- filter(mt, species %in% c("gp17.5"))
  mt_lysin <- filter(mt, species %in% c("gp3.5"))
  
  
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
    gp16 <- mt %>% filter (species %in% c("gp16"), protein>=19744)
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




