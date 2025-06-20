#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=10G
#SBATCH -J check_raws_cluster_name
#SBATCH -o slurm_outputs/check_raws_cluster_name-%j.out
#SBATCH -e slurm_outputs/check_raws_cluster_name-%j.err
#SBATCH --array=0-99

# I've done a very sloppy parallelization of this, use a jobarray and from there let python decide
# whether to open the file (since that makes up almost all of the runtime) or 
# to skip it because its another tasks problem

# defining variables
CLN="cluster_name"
LOAD_LSST="load_pipeline_path" # REPLACE WITH load_pipeline_path FOR TEMPLATE
CLUSTER_DIR="cluster_dir"
PY_SCRIPTS="py_scripts"

# navigate to cluster directory
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# in-case a previous list/out.txt exists, delete it
# this may cause problems later if this script is used multiple times... watch out for that
rm -r raws/corrupt_raws_${SLURM_ARRAY_TASK_ID}.txt
rm -r raws/out_${SLURM_ARRAY_TASK_ID}.txt
mkdir -p corrupt_raws

# begin checking each raw file
for file in raws/*.fz; do
	echo "Verifying ${file}"
	python python_scripts/check_data.py "${file}" ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_MAX} > raws/out_${SLURM_ARRAY_TASK_ID}.txt
	# not-empty output imples an intact file!
	if [ -s raws/out_${SLURM_ARRAY_TASK_ID}.txt ]; then
		rm -r raws/out_${SLURM_ARRAY_TASK_ID}.txt
		echo "${file} is intact!"
		continue
	else
		# out.txt is empty, implies file is corrupt so separate it from the rest!
		echo "${file} is corrupt! Writing to corrupt_raws ..."
		echo "${file}" >> raws/corrupt_raws_${SLURM_ARRAY_TASK_ID}.txt
		# mv "${file}" corrupt_raws/
		# write a script to move files in this list elsewhere
	fi
done
