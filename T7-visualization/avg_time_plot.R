library(plyr)
library(tidyverse)
library(ggplot2)

#get filenames
filenames = list.files(pattern = "*.tsv")

# load all files into a single data frame and round timepoints to nearest 5. Group data by time and species
cts = ldply(filenames, read_tsv) %>% group_by(time, species)
cts$time = round(cts$time/5)*5

#set proteins of interest
species.of.interest = c("gp0.7", "gp1", "gp10A", "gp3.5")

#filter for species of interest and find average amount of protein for each species at each timepoint
cts.to.graph = cts %>% filter(species %in% species.of.interest) %>% dplyr::summarise(avg = mean(protein), time = time, species = species)

#graph
time_plot = ggplot(cts.to.graph, aes_string(x="time", y="avg", color="species")) + 
  geom_line() +
  theme_bw()
time_plot

#make nice html display version of graph
library(plotly)
pPlotly <- ggplotly(time_plot)
htmlwidgets::saveWidget(pPlotly, "time_plot.html")