library(plyr)
library(tidyverse)
library(ggplot2)

#**********************************************************************
#BEFORE RUNNING ENTER THE FOLLOWING THINGS
# wt_wd = "[Enter wildtype working directory file path]"
# mt_name = "[Enter name of mutant]"
# mt_wd = "[Enter mutant working directory file path"]
#**********************************************************************
#TO RESTRICT GRAPHS TO SPECIFIC GENE CLASS ENTER THE FOLLOWING THINGS
# class.number = [1, 2, or 3]
#**********************************************************************

#set wt
setwd(wt_wd)

#get file names
filenames = list.files(pattern = "*.tsv")

#load all files into a single data frame and round timepoints to nearest 5. Group data by species
wt.cts = ldply(filenames, read_tsv) %>% group_by(species)
wt.cts$time = round(wt.cts$time/5)*5

#load gene display file
disp = read.csv("gene_display.csv", stringsAsFactors = FALSE)

wt.cts$protein[wt.cts$species == "gp1"] = wt.cts$protein[wt.cts$species == "gp1"] + wt.cts$protein[wt.cts$species == "gp1+gp3.5"]


#find average amount of protein for each species at final timepoint
#input 
#class.number = 1, 2, or 3
#followed by 
#class.number == 1, 2, or 3
#to specify boxplot to only display gps of one class
wt.cts.to.graph = wt.cts %>% filter(time == max(wt.cts$time))

wt.cts.to.graph = wt.cts.to.graph %>% distinct() %>% left_join(disp, by = "species")

disp = disp %>% arrange(sort.order) %>% subset.data.frame(gene.class == class.number)

wt.cts.to.graph$species = factor(wt.cts.to.graph$species, levels = disp$species)

#omit NA
wt.cts.to.graph <- na.omit(wt.cts.to.graph)

#set mutant wd
setwd(mt_wd)

#get file names
filenames = list.files(pattern = "*.tsv")

#load all files into a single data frame and round timepoints to nearest 5. Group data by species
cts = ldply(filenames, read_tsv) %>% group_by(species)
cts$time = round(cts$time/5)*5

#load gene display file
disp = read.csv("gene_display.csv", stringsAsFactors = FALSE)

cts$protein[cts$species == "gp1"] = cts$protein[cts$species == "gp1"] + cts$protein[cts$species == "gp1+gp3.5"]


#find average amount of protein for each species at final timepoint
cts.to.graph = cts %>% filter(time == max(cts$time))

cts.to.graph = cts.to.graph %>% distinct() %>% left_join(disp, by = "species")

disp = disp %>% arrange(sort.order) %>% subset.data.frame(gene.class == class.number)

cts.to.graph$species = factor(cts.to.graph$species, levels = disp$species)

#test for significance
sig = list()
species_gp = c(disp$species)
for (i in seq(1:52)) {
  wt.species = wt.cts.to.graph %>% filter(species == disp$species[i])
  mt.species = cts.to.graph %>% filter(species == disp$species[i])
  sig = sig %>% append(p.adjust(wilcox.test(wt.species$protein, mt.species$protein)$p.value))
  significant = data.frame(species= species_gp, p_values = cbind(sig[1:52]), color = "red")
}

for (i in seq(1:52)) {
  if (sig[i] > 0.05){significant$color[i] = "#A9A9A9"}
}

cts.to.graph = na.omit(cts.to.graph)
cts.to.graph$color = "red"
for (i in seq (1:nrow(cts.to.graph))){
  for (j in seq(1:52)) {
    if (cts.to.graph$species[i] == significant$species[j])
      cts.to.graph$color[i] = significant$color[j]
  }
}

significant$wtcolor = ""
for (i in seq(1:52)) {
  if (significant$color[i] == "#A9A9A9"){
    significant$wtcolor[i] = "black"
  }
  if (significant$color[i] == "red"){
    significant$wtcolor[i] = "blue"
  }
}

wt.cts.to.graph$color = ""
for (i in seq (1:nrow(wt.cts.to.graph))){
  for (j in seq(1:52)) {
    if (wt.cts.to.graph$species[i] == significant$species[j])
      wt.cts.to.graph$color[i] = significant$wtcolor[j]
  }
}

#new boxplot coloring
mutantcompare_genomes_boxplot <- ggplot(na.omit(cts.to.graph), aes(x= species, y = protein), legend(x = 250, y = 200, legend = c("Wildtype", mt_name), fill = c("red", "blue"))) + 
  stat_boxplot( geom ='errorbar', aes(color = cts.to.graph$color), alpha = 0.5, width = 0.7) +
  geom_boxplot(fill = significant$color, aes(color = cts.to.graph$color), alpha = 0.5, outlier.size= 0.5, width = 0.7)  +
  geom_boxplot(data = wt.cts.to.graph,fill = significant$wtcolor, aes(color = wt.cts.to.graph$color), alpha = 0.5, outlier.size = 0.5, width = 0.4)+
  stat_boxplot( data = wt.cts.to.graph, geom ='errorbar', aes(color = wt.cts.to.graph$color), alpha =0.5, width = 0.4) +
  scale_y_continuous(trans='log10') +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_identity(name = "Legend",
                       breaks = c("blue", "red"),
                       labels = c("Wildtype", "3.5 Mutant"),
                       guide = "legend")
mutantcompare_genomes_boxplot
