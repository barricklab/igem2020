#Instructions for running T7 simulations on TACC (ls5)
# 2020-08-04

##########################################################################
# Step 1: Set up paths and modules needed on TACC in our .bashrc

# Edit ~/.bashrc
nano ~/.bashrc

# Add these lines under "PLACE MODULE COMMANDS HERE and ONLY HERE."
module load git
module load Rstats
module load python3

# Add these lines under "PLACE Environment Variables including PATH here"
export PATH="$PATH:$WORK/src/igem2020/T7-simulation"
export PS1='LS5:\w\$ '

# Now log in again or `source $HOME/.bashrc` so that the changes take effect.

##########################################################################
# Step 2: Install PineTree and our code

# Create directory where we will store all code
mkdir -p $WORK/src

# Clone PineTree repo to install latest version
cd $WORK/src
git clone https://github.com/clauswilke/pinetree.git
cd pinetree
python3 setup.py install --user
# Note: the pip3 install command on GitHub doesn't work on ls5 for installation

# Clone iGEM 2019 repo
cd $WORK/src
git clone https://github.com/barricklab/igem2020.git


##########################################################################
# Step 3: Run simulations

# Create a directory containing the genome of interest and job control files

# Create directory where we will store all code
mkdir -p $SCRATCH/T7-example-run

# Place the genome file there (we will copy the wild-type as an example)
cp $WORK/src/igem2020/T7-simulation/T7-WT/T7_genome.gb $SCRATCH/T7-example-run

# Copy the job script that is checked into the repository
cp $WORK/src/igem2020/T7-simulation/TACC/launcher.slurm $SCRATCH/T7-example-run

#Change into the job directory
cd $SCRATCH/T7-example-run
