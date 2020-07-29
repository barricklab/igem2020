require(ggplot2)
require(tidyverse)
require(plotly)
#read in both the wild type and mutant phage count files, files names may vary 
wt = read_tsv("T7-WT-3255.counts.tsv")
mt = read_tsv("T7_replace_10_with_sfGFP.counts.tsv")
disp = read.csv("gene_display.csv", stringsAsFactors = FALSE)

#filter for all of the genes at the maximum time point of the simulation 
wt.to.graph = wt %>% filter(species %in% disp$species, time == max(wt$time)) %>% left_join(disp, by = "species")

mt.to.graph = mt %>% filter(species %in% disp$species, time == max(mt$time)) %>% left_join(disp, by = "species")

disp = disp %>% arrange(sort.order)

wt.to.graph$species = factor(wt.to.graph$species, levels = disp$species) 
wt.to.graph = wt.to.graph %>% complete(species)

mt.to.graph$species = factor(mt.to.graph$species, levels = disp$species)
mt.to.graph = mt.to.graph %>% complete(species)


#creates an overlapping bar plot for the mutant and wild-type genes, graphing maximum protein production of each  
#you can change the y axis to other variables such as "transcript", or "ribo-density"

compare_genomes_barplot <- ggplot(NULL, aes(x = species, y = protein)) + 
  geom_bar(aes(fill = "Mutant"), data = mt.to.graph, alpha = 0.5, stat = "identity") +
  geom_bar(aes(fill = "WT"), data = wt.to.graph, alpha = 0.5, stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90)) 
compare_genomes_barplot

#generates plotly for the graph and saves it locally in the working directory 
library(plotly)
pPlotly <- ggplotly(compare_genomes_barplot)
htmlwidgets::saveWidget(pPlotly, "compare_genomes.html")
