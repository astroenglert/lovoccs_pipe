#!/bin/bash

# The commands for slurm start with #SBATCH
# All slurm commands need to come before the program 
# you want to run.  In this example, 'echo "Hello World!"
# is the command we are running.
#
# This is a bash script, so any line that starts with # is
# a comment.  If you need to comment out an #SBATCH line 
# use ##SBATCH 
#
# To submit this script to slurm do:
#    sbatch batch.script
#
# Once the job starts you will see a file DownloadNOAO-****.out
# The **** will be the slurm JobID

# --- Start of slurm commands -----------

# Request an hour of runtime:
#SBATCH --time=23:59:59

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH --mem=15G

# Specify a job name:
#SBATCH -J initialize_repo_cluster_name

# Specify an output file
# %j is a special variable that is replaced by the JobID when 
# job starts
#SBATCH -o slurm_outputs/initialize_repo_cluster_name-%j.out
#SBATCH -e slurm_outputs/initialize_repo_cluster_name-%j.out

#----- End of slurm commands ----

# source LSST_LOC/loadLSST.bash
setup lsst_distrib

currentDir=${PWD}

# creating butler_repo (if it doesn't already exist)
mkdir -p "${currentDir}/butler_repo"

# move into the butler_repo folder to make the next steps a little easier
cd "${currentDir}/butler_repo"

# tell LSP to create a repository here and register it for DECam data
butler create butler_repo
butler register-instrument butler_repo lsst.obs.decam.DarkEnergyCamera

# copy a handful of calibration frames from LSP to this repository
# unfortunately no soft-link option for this, but these calibration frames aren't that big
butler write-curated-calibrations butler_repo lsst.obs.decam.DarkEnergyCamera







