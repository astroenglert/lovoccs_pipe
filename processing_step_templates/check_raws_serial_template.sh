#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=15G
#SBATCH -J check_raws_cluster_name
#SBATCH -o slurm_outputs/check_raws_cluster_name-%j.out
#SBATCH -e slurm_outputs/check_raws_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir"
PY_SCRIPTS="py_scripts"

# navigate to cluster directory
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# in-case a previous list exists, delete it
# this may cause problems later if this script is used multiple times... watch out for that
rm -r raws/corrupt_raws.txt
mkdir -p corrupt_raws

# begin checking each raw file
for file in raws/*.fz; do
	echo "Verifying ${file}"
	python python_scripts/check_data.py "${file}" > raws/out.txt
	if [ -s raws/out.txt ]; then
		# not empty output, implies intact file, so continue to next file
		rm -r raws/out.txt
		echo "${file} is intact!"
		continue
	else
		# out.txt is empty, implies file is corrupt so separate it from the rest!
		echo "${file} is corrupt! Moving to corrupt_raws ..."
		echo "${file}" >> corrupt_raws/corrupt_raws.txt
		mv "${file}" corrupt_raws/

		# update raw_info csv
		echo "${file:4}"
		sed "*${file:4}*/s/FALSE/TRUE/" raws/raw_info.csv > test_update.csv
	fi
done






