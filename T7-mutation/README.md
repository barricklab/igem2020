#  Setup

Install _breseq_ using one of the methods described in the [_breseq_ instructions](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/installation.html)

**Note:** Because we are only using the `gdtools` utility command from _breseq_. You do not need to install the requirements for _breseq_ (bowtie2 and R)!

# Create GenomeDiff file of mutations

Create a GenomDiff file describing the mutations you want to 'apply' to the genome.


This one `delete_1.4-1.6.gd` deletes genes 0.4, 1.5, 1.6:
```txt
#=GENOME_DIFF	1.0
DEL	.	.	T7	7608	559
```

This one `replace_0.5_with_sfGFP.gd` deletes gene 0.5 and inserts sfGFP in its place:
```txt
#=GENOME_DIFF	1.0
INT	.	.	T7	1278	156	sfGFP:1-717
DEL	.	.	sfGFP	1	717
```

In these `*.gd` files, it is important that the columns are separated by 'tabs' and not spaces!

For a full description of the syntax, including other types of mutations you can use, see the [GenomeDiff Format reference](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_format.html).

# Generate mutated genome GenBank file

Use the `gdtools` utility from _breseq_ to mutate the sequence:
```sh
gdtools APPLY -f GENBANK -r T7_genome.gb -o T7_delete_1.4-1.6.gb delete_1.4-1.6.gd
gdtools APPLY -f GENBANK -r T7_genome.gb -r sfGFP.gb -o T7_replace_0.4_with_sfGFP.gb replace_0.4_with_sfGFP.gd
```

This creates the new GenBank files `T7_delete_1.4-1.6.gb` and `T7_replace_0.4_with_sfGFP.gb`. Take a look at them to see if they've been mutated!

If you want to understand the options here, type just `gdtools APPLY` at the command line to display the help.

# More information

* https://github.com/barricklab/breseq
