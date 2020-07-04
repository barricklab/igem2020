#  Setup

Create a virtual environment and install requirements
```sh
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

When you start in the future, you will only need to run this command:
```
source env/bin/activate
```

# Running the wild-type model

```sh
cd T7-WT
python3 phage_model.py (genbank input filepath) [output folder]
```
This will take more than ten minutes to run! Output is in the file `phage_counts.tsv`.
You can examine progress and the results while the simulation is still running,
for example using:

```sh
tail phage_counts.tsv (genbank input filepath) [output folder]
```

If you are impatient, you can examine the example TSV output file that is in the repository and use it to visualize the results in the next step: [`T7-WT/example_output_phage_counts.tsv`](T7-WT/example_output_phage_counts.tsv).

# Visualizing output

PineTree output is a data table in [tidy format](https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html).

You can filter rows corresponding to different genes/species and graph them. For example, these commands in R will graph the concentration of a few genes over time:

```R
require(ggplot2)
require(tidyverse)

cts = read_tsv("example_output_phage_counts.tsv")

species.of.interest = c("gene 0.3", "rnapol-1", "rnapol-3.5", "gene 6", "gene 9")
cts.to.graph = cts %>% filter(species %in% species.of.interest)

#Choose from protein, transcript, or ribo_density
to.graph="protein"

#Create plot
ggplot(cts.to.graph, aes_string(x="time", y=to.graph, color="species")) +
  geom_line() +
  theme_bw()
```

# More information

* https://pinetree.readthedocs.io/en/latest/intro.html
* https://github.com/alexismhill3/pinetree-demo
* https://github.com/clauswilke/pinetree
