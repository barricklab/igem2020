#  Setup

**You must be using a newer version of _breseq_ (v0.35.3rc1 or later) for these instructions to work!**

Use this link to find this version: https://github.com/barricklab/breseq/releases/tag/v0.35.3rc1

Install _breseq_ using one of the methods described in the [_breseq_ instructions](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/installation.html)

**Note:** Because we are only using the `gdtools` utility command from _breseq_. You do not need to install the requirements for _breseq_ (bowtie2 and R)!

# Create GenomeDiff file of mutations

Create a GenomDiff file describing the mutations you want to `APPLY` to the genome.


This one `delete_1.4-1.6.gd` deletes genes 0.4, 1.5, 1.6:
```txt
#=GENOME_DIFF	1.0
DEL	.	.	T7	7608	559
```
**Note:** This `DEL` removes 559 bases starting at position 7608. That is, it deletes bases 7608-8166.

This one `replace_0.4_with_sfGFP.gd` deletes gene 0.4 and inserts sfGFP in its place:
```txt
#=GENOME_DIFF	1.0
INT	.	.	T7	1278	156	sfGFP:1-717
```
**Note:** This `INT` replaces 156 bases starting at position 1278 with the specified region from the other sequence specified as `seq_id:start-end`. This means that bases 1278-1433 are replaced in this example.

This one `insert_sfGFP_after_0.3.gd` inserts sfGFP between genes 0.3 and 0.4:
```txt
#=GENOME_DIFF	1.0
DEL	1	.	T7	925	509
INT	2	.	T7	924	0	gp0.3:1-354
INT	3	.	T7	924	0	sfGFP:1-717
INT	4	.	T7	924	0	gp0.4:1-156
```
**Note:** In this case, the genes overlap by a base in the original genome, so it is necessary to first delete them, and then add back 0.3, sfGFP, and 0.4. an `INT` with a size of zero will insert the new region after the given position.

In all of these these `*.gd` files, it is important that the columns are separated by 'tabs' and not spaces!

For a full description of the syntax, including other types of mutations you can use, see the [GenomeDiff Format reference](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_format.html).

# Generate mutated genome GenBank file

Use the `gdtools` utility from _breseq_ to create the mutated T7 genomes using these input files:
```sh
gdtools APPLY -s T7 -f GENBANK -r T7_genome.gb -o T7_delete_1.4-1.6.gb delete_1.4-1.6.gd
gdtools APPLY -s T7 -f GENBANK -r T7_genome.gb -r sfGFP.gb -o T7_replace_0.4_with_sfGFP.gb replace_0.4_with_sfGFP.gd
gdtools APPLY -s T7 -f GENBANK -r T7_genome.gb -r sfGFP.gb -r T7_genes.gb -o T7_insert_sfGFP_after_0.3.gb insert_sfGFP_after_0.3.gd
```

This creates the new GenBank files `T7_delete_1.4-1.6.gb`, `T7_replace_0.4_with_sfGFP.gb`, `T7_insert_sfGFP_after_0.3.gb`. Take a look at them to see if they've been mutated!

**Notes:** There are additional `-r` entries when a file is needed to provide the sequence and annotation for an `INT`. The `-s` option makes sure that only the mutated T& is included in the output, rather than also including the sfGFP sequence, for example.

If you want to further understand the options here, type just `gdtools APPLY` at the command line (adding no other options) to display the help.

# More information

* https://github.com/barricklab/breseq
