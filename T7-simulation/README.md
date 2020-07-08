#  Setup

Create a virtual environment and install required Python modules:
```sh
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

When you start in the future, you will only need to run this command:
```
source env/bin/activate
```

# Running the T7 phage model

If you run the script `phage_model.py` with the `--help` option, then it will print the help that describes all of the command line arguments.

```sh
python3 phage_model.py  --help
```

Here is an example of how to run the wild-type phage model for a simulation lasting 900 seconds (15 minutes) after a phage begins infecting an _E. coli_ cell.

**Important:** The following commands are written to be run from a terminal window in which your current working directory is the main `T7-simulation` directory.

```sh
python3 phage_model.py -t 900 -i T7-WT/T7_genome.gb -o user-output/T7-WT
```
This will take more than five minutes to run! Simulation output is in the file `user-output/T7-WT.counts.tsv` and a log file showing how the genome was loaded and what version of the script you used is in the file `user-output/T7-WT.log`. You can open either of these in a text editor.

You can examine progress and the results while the simulation is still running, for example by opening a new terminal window and using this command to continuously show you the last lines of the output file:

```sh
tail -f user-output/T7-WT.counts.tsv
```

**Tip:** You can press `control-C` in a terminal window to kill either of these commands and get back to the prompt if you decide you don't want to wait.

Also, if you are impatient, you can examine the example TSV output file that is in the repository and use it to visualize the results in the next step: [`T7-simulation/T7-WT/T7.example.counts.tsv`](T7-simulation/T7-WT/T7.example.counts.tsv).

Mutated phage genomes that you can also test simulating are provided in the `T7-mutant` directory.

You can store new mutant genomes that you create for testing in `user-input` and send the output to `user-output` so that these won't be tracked on GitHub.

# Visualizing output

PineTree output is a data table in [tidy format](https://cran.r-project.org/web/packages/tidyr/vignettes/tidy-data.html).

You can filter rows corresponding to different genes/species and graph them. For example, these commands in R will graph the concentration of a few genes over time:

```R
require(ggplot2)
require(tidyverse)

cts = read_tsv("T7.example.counts.tsv")

species.of.interest = c("gene 0.3", "rnapol-1", "rnapol-3.5", "gene 6", "gene 9")
cts.to.graph = cts %>% filter(species %in% species.of.interest)

#Choose from protein, transcript, or ribo_density
to.graph="protein"

#Create plot
ggplot(cts.to.graph, aes_string(x="time", y=to.graph, color="species")) +
  geom_line() +
  theme_bw()
```

R scripts for visualization and analysis are provided in [`T7-visualization`](T7-visualization).

# More information

* https://pinetree.readthedocs.io/en/latest/intro.html
* https://github.com/alexismhill3/pinetree-demo
* https://github.com/clauswilke/pinetree
