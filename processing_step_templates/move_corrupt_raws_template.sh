#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=10G
#SBATCH -J move_corrupt_raws_cluster_name
#SBATCH -o slurm_outputs/move_corrupt_raws_cluster_name-%j.out
#SBATCH -e slurm_outputs/move_corrupt_raws_cluster_name-%j.err

# defining variables
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir"
PY_SCRIPTS="py_scripts"

# navigate to cluster directory
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
# we don't even need to for this step, since its just shuffling files around

#source ${LOAD_LSST}
#setup lsst_distrib

# begin checking each raw file
# this little trick won't work if the files have whitespace in their names
# but I don't think that ever happens so it shouldn't cause any issues
for file in raws/corrupt_raws_*; do
	echo "Moving raws in ${file}"
	while read line; do
		mv ${line} corrupt_raws/$(basename ${line})
	done<${file}
	mv ${file} corrupt_raws/$(basename ${file})
done

