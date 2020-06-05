require(ggplot2)
require(tidyverse)

cts = read_tsv("phage_counts.tsv")

species.of.interest = c("gene 0.3", "rnapol-1", "rnapol-3.5", "gene 6", "gene 9")
cts.to.graph = cts %>% filter(species %in% species.of.interest)

#Choose from protein, transcript, or ribo_density
to.graph="protein"

#Create plot
ggplot(cts.to.graph, aes_string(x="time", y=to.graph, color="species")) + 
  geom_line() +
  theme_bw()

       