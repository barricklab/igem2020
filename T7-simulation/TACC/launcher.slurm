#!/bin/bash
#
# Modified by Jeff Barrick for T7 simulations
#
# Laucher batch script file for TACC systems (like Frontera, Stampede2, etc.)
# Si Liu
# July 13, 2020
#
# Simple SLURM script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the LAUNCHER_JOB_FILE to specify 
# each executable per process.
#-------------------------------------------------------
# 
#
#
# Option -t time is HH::MM:SS
# Option -p should be normal|development
# Option should be 48 x N
# Option should be the the number of commands/48 rounded UP to an integer.
#
#         <------ Setup Parameters ------>
#
#SBATCH -J T7-pinetree
#SBATCH -A breseq
#SBATCH -N 1                # Must be commands/48 rounded UP to an integer
#SBATCH -n 48               # Option should be 48 x N
#SBATCH -p development      # Use normal if N>48 or -t > 1 hr
#SBATCH -o Parametric.%j.out
#SBATCH -e Parametric.%j.err
#SBATCH -t 01:00:00         # Format HH::MM:SS
#SBATCH --mail-user=username@tacc.utexas.edu
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

#------------------------------------------------------

module load launcher
export LAUNCHER_WORKDIR=.
export LAUNCHER_JOB_FILE=commands.sh

$LAUNCHER_DIR/paramrun

