library(plyr)
library(tidyverse)
library(ggplot2)
source("bs_lt_functions.R")

#get filenames
filenames = list.files(pattern = "*.tsv")

#read gene display file
disp = read.csv("gene_display.csv", stringsAsFactors = FALSE)

#create df to hold the proteins of each gene of each trial at lysis time of trial
protein_counts= data.frame()  

for (i in filenames){ #find lysis time of each trial, add protein counts of each gene to protein_counts df
  file_name = i
  mt = read_tsv(file_name)
  lt = find_lts(mt)
  mt$protein[mt$species == "gp1"] = mt$protein[mt$species == "gp1"] + mt$protein[mt$species == "gp1+gp3.5"]
  if (lt != -1){ 
    prot <- filter(mt, species %in% disp$species, time == lt)
    protein_counts = rbind(protein_counts, prot)}
} 

#create vector for gene name (s) and for average protein of that gene (p) -> will be converted into df for plotting later
s <- c()
p <- c()
g <- disp$gene.class
o <- disp$sort.order
l <- c()
u <- c()
for (ss in disp$species){ #for every gene listed in gene display file, filter for gene name, find average protein count
                         #add gene name to s and average protein count to p
  
  a <- protein_counts %>% filter(species == ss)
  avg_sp = mean(a$protein)
  lq <- quantile(a$protein, probs = c(0.025))
  uq <- quantile(a$protein, probs = c(0.975))
  l <- append (l, lq)
  u <- append (u, uq)
  s <- append(s, ss)
  p <- append(p, avg_sp)
  
}
#create df that takes in s and p vector for plotting
avg_pc = data.frame(s, p, o, g, l, u)
colnames(avg_pc) = c("species", "average", "sort.order", "gene.class", "lower95", "upper95")


disp = disp %>% arrange(sort.order)
avg_pc$species = factor(avg_pc$species, levels = disp$species)

#plots 
compare_genomes_barplot <- ggplot(avg_pc, aes(x = species, y = average, fill = as.factor(gene.class))) + 
  geom_bar( alpha = 0.5, stat = "identity") +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("#CC79A7", "#56B4E9", "#E69F00"), name = "Class", labels = c("I", "II", "III")) +
  theme(axis.text.x = element_text(angle = 90), axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank())  
compare_genomes_barplot


#generates plotly for the graph and saves it locally in the working directory 
library(plotly)
pPlotly <- ggplotly(compare_genomes_barplot)

htmlwidgets::saveWidget(pPlotly, "Holin_replace_8-AE.html")