#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=15G
#SBATCH -J ingest_data_cluster_name
#SBATCH -o slurm_outputs/ingest_data_cluster_name-%j.out
#SBATCH -e slurm_outputs/ingest_data_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name" 
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # REPLACE WITH CLUSTER DIRECTORY EVENTUALLY
PY_SCRIPTS="py_scripts"
#TODO Pass this via run_steps rather than hard-coded
CALIB_CATALOG_REPO="/gpfs/data/idellant/Clusters/calib_catalog_repo" 

# navigate to cluster directory
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# import reference calibs, using symlink to avoid fussing over calibration dates

echo "Importing calibration frames!"
butler import --transfer direct repo/repo "${CALIB_CATALOG_REPO}"/import_ready_calibs_7

#2024-09 AE: Consolidated calibrations into a single collection, no need to import these separately anymore

# butler import --transfer direct repo/repo "${CALIB_CATALOG_REPO}"/import_ready_skyframes
# butler import --transfer direct repo/repo "${CALIB_CATALOG_REPO}"/BF-Effect/import_ready_bfk_2

for catalog in "${CALIB_CATALOG_REPO}/gen3_formatted_new/${CLN}/"*; do

	# basename gets the title of the catalog without any of the filepath
	# BASE_NAME=$(basename ${catalog%.*})
	BASE_NAME=$(basename ${catalog})
	
	echo "Importing ${BASE_NAME}!"
	# create a dataset-type for the catalog
	butler register-dataset-type repo/repo ${BASE_NAME} SimpleCatalog htm7

	# now ingest the formatted catalog into the refcats collection
	# to reference this catalog point to the refcats collection and specify the dataset_type... 
	# the dataset_type is effectively the title of each folder in .../calib_catalog_repo/gen3_formatted/cluster_name/*
	butler ingest-files -t symlink repo/repo ${BASE_NAME} refcats ${CALIB_CATALOG_REPO}/gen3_formatted_new/${CLN}/${BASE_NAME}/filename_to_htm.ecsv

done

echo "Importing raw frames!"

# ingest the raw data
butler ingest-raws repo/repo raws/*.fits.fz --transfer direct

echo "Creating visit info!"

butler define-visits repo/repo lsst.obs.decam.DarkEnergyCamera


