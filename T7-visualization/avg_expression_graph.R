library(plyr)
library(tidyverse)
library(ggplot2)

#get filenames
filenames = list.files(pattern = "*.tsv")

# load all files into a single data frame and round timepoints to nearest 5. Group data by time and species
cts = ldply(filenames, read_tsv) %>% group_by(species)
cts$time = round(cts$time/5)*5

# load gene display file
disp = read.csv("gene_display.csv", stringsAsFactors = FALSE)

#find average amount of protein for each species at final timepoint
cts.to.graph = cts %>% filter(species %in% disp$species, time == max(cts$time)) %>% dplyr::summarise(avg = mean(protein), species = species, lower95 = quantile(protein, probs = c(0.025)), upper95 = quantile(protein, probs = c(0.975)))

cts.to.graph = cts.to.graph %>% distinct() %>% left_join(disp, by = "species")

disp = disp %>% arrange(sort.order)

cts.to.graph$species = factor(cts.to.graph$species, levels = disp$species) 

#plots
compare_genomes_barplot <- ggplot(cts.to.graph, aes(x = species, y = avg)) + 
  geom_bar(aes(fill = "WT"), alpha = 0.5, stat = "identity") +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.5, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) 
compare_genomes_barplot

#generates plotly for the graph and saves it locally in the working directory 
library(plotly)
pPlotly <- ggplotly(compare_genomes_barplot)

htmlwidgets::saveWidget(pPlotly, "WT-profile.html")