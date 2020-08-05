#  T7 Data Visualization

For scripts that visualize output data from the T7 simulation

The `time_course.R` script plots the abundance of selected proteins over time from one simulation.

The `mutation_comp_ordered.R` script generates an overlapping barplot that compares the protein production levels of all genes from two genomes of interest. The input file is set by creating the variable  "input.prefix" in the RStudio console and setting it equal to the name of your input file minus the .counts.tsv suffix. For example, if you want to visualize the file `T7_delete_10.counts.tsv` you would run the command `input.prefix = "T7_delete_10"` and then run the entire script. Currently the wild type file for comparison is hardcoded to `T7-WT_24.counts.tsv` so make sure that, `gene_display.csv`, and your mutated count file are all in the working directory.

