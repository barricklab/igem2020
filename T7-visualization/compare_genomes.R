require(ggplot2)
require(tidyverse)
require(plotly)
#read in both the wild type and mutant phage count files, files names may vary 
wt = read_tsv("example_output_phage_counts.tsv")
mt = read_tsv("mutant_phage_counts.tsv")

#filter for all of the genes at the maximum time point of the simulation 
wt.to.graph = wt %>% filter((str_detect(wt$species, "^gene")) , time == max(wt$time))
mt.to.graph = mt %>% filter((str_detect(mt$species, "^gene")), time == max(mt$time))

#creates an overlapping bar plot for the mutant and wild-type genes, graphing maximum protein production of each  
#you can change the y axis to other variables such as "transcript", or "ribo-density"
compare_genomes_barplot <- ggplot(NULL, aes_string(x="species", y = "protein")) + 
  geom_bar(aes(fill = "Mutant"), data = mt.to.graph, alpha = 0.5, stat = "identity") +
  geom_bar(aes(fill = "WT"), data = wt.to.graph, alpha = 0.5, stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90)) 
compare_genomes_barplot

#generates plotly for the graph and saves it locally in the working directory 
library(plotly)
pPlotly <- ggplotly(compare_genomes_barplot)
htmlwidgets::saveWidget(pPlotly, "compare_genomes.html")

