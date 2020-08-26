library(plyr)
library(tidyverse)
library(ggplot2)

# Make sure to specify output.data as "protein" or "transcript" BEFORE you run this script! 

# get filenames
filenames = list.files(pattern = "*.tsv")

# load all files into a single data frame and round timepoints to nearest 5. Group data by time and species
cts = ldply(filenames, read_tsv) %>% group_by(species)
cts$time = round(cts$time/5)*5

# load gene display file
disp = read.csv("gene_display.csv", stringsAsFactors = FALSE)

if (output.data == "protein"){
  #fix gp1 count
  cts$protein[cts$species == "gp1"] = cts$protein[cts$species =="gp1"] + cts$protein[cts$species == "gp1+gp3.5"]

  #find average amount of protein for each species at final timepoint
  cts.to.graph = cts %>% filter(species %in% disp$species, time == max(cts$time)) %>% dplyr::summarise(avg = mean(protein), species = species, lower95 = quantile(protein, probs = c(0.025)), upper95 = quantile(protein, probs = c(0.975)))

  cts.to.graph = cts.to.graph %>% distinct() %>% left_join(disp, by = "species")
}

if (output.data == "transcript"){

  #find average amount of protein for each species at final timepoint
  cts.to.graph = cts %>% filter(species %in% disp$species, time == max(cts$time)) %>% dplyr::summarise(avg = mean(transcript), species = species, lower95 = quantile(transcript, probs = c(0.025)), upper95 = quantile(transcript, probs = c(0.975)))
  
  cts.to.graph = cts.to.graph %>% distinct() %>% left_join(disp, by = "species")
}



disp = disp %>% arrange(sort.order)

cts.to.graph$species = factor(cts.to.graph$species, levels = disp$species) 

#plots

compare_genomes_barplot <- ggplot(cts.to.graph, aes(x = species, y = avg, fill = as.factor(gene.class))) + 
  geom_bar(alpha = 0.5, stat = "identity") +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#E69F00"), name = "Class", labels = c("I", "II", "III")) +
  theme(axis.text.x = element_text(angle = 90), axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
compare_genomes_barplot

#generates plotly for the graph and saves it locally in the working directory
#make sure to rename plotly file below
library(plotly)
pPlotly <- ggplotly(compare_genomes_barplot)

htmlwidgets::saveWidget(pPlotly, "MT-3.5-profile+RNAse-transcript.html")
