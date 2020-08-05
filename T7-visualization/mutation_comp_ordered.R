require(ggplot2)
require(tidyverse)
require(plotly)

##shell stuff##
#!/usr/bin/env Rscript
#add the name of the mutant file, without the .counts.tsv extension part
#ex: input.prefix = "T7_replace_10_with_sfGFP"
#after you define the line above, you should be able to run the script in one swift go, assuming you have the wt and disp files in your directory 
option.input.prefix = input.prefix
mt = data.frame()
######################################
# check that file exists 
if (file.exists(paste0(option.input.prefix, ".counts.tsv"))) {
  mt <- read_tsv(paste0(option.input.prefix, ".counts.tsv"), col_names=TRUE )
} 


# Default to the same output prefix as input prefix. 
if (exists("output.prefix")) {
  option.output.prefix = output.prefix
} else {
  option.output.prefix = input.prefix
}

################################


#read in both the wild type and mutant phage count files, files names may vary 
wt = read_tsv("T7-WT_24.counts.tsv")

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

htmlwidgets::saveWidget(pPlotly, "delete_10_compare_genomes.html")

htmlwidgets::saveWidget(pPlotly, "compare_genomes.html")
