#!/bin/bash

# This is an example batch script for slurm on Oscar
# 
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
#SBATCH -J download_noao_cluster_name

# Specify an output file
# %j is a special variable that is replaced by the JobID when 
# job starts
#SBATCH -o slurm_outputs/download_noao_cluster_name-%j.out
#SBATCH -e slurm_outputs/download_noao_cluster_name-%j.out

#----- End of slurm commands ----

# Run a command

download_command_filename="processing_step/download_raw_cluster_name.sh"

noao_query_py="template_folder/noao_query_raw.py"
raw_image_folder="raw"
current_path=${PWD}


#-----------------------
module load curl 

source ${HOME}/.bashrc.ext

#[ ! -d /gpfs/data/idellant/Clusters/cluster_name ] && mkdir /gpfs/data/idellant/Clusters/cluster_name
#cd /gpfs/data/idellant/Clusters/cluster_name

[ -f "${download_command_filename}" ] && rm "${download_command_filename}"

python ${noao_query_py} cluster_name ${download_command_filename}
bash ${download_command_filename}

[ ! -d ${raw_image_folder} ] && mkdir ${raw_image_folder} 
mv *.fz ${raw_image_folder}/


#cd ${current_path}
deactivate
