#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=15G
#SBATCH -J gotta_blast_cluster_name
#SBATCH -o slurm_outputs/gotta_blast_cluster_name-%j.out
#SBATCH -e slurm_outputs/gotta_blast_cluster_name-%j.err

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

# the following datasets are removed from the repository
echo Blasting postISRCCD
butler prune-datasets repo/repo --find-all --unstore --datasets postISRCCD --no-confirm DECam/processing/*
echo Blasting icSrc
butler prune-datasets repo/repo --find-all --unstore --datasets icSrc --no-confirm DECam/processing/*
echo Blasting icExp
butler prune-datasets repo/repo --find-all --unstore --datasets icExp --no-confirm DECam/processing/*
echo Blasting calexps
butler prune-datasets repo/repo --find-all --unstore --datasets calexp --no-confirm DECam/processing/*
echo Blasting icExpBackground
butler prune-datasets repo/repo --find-all --unstore --datasets icExpBackground --no-confirm DECam/processing/*
echo Blasting icSrcSchema
butler prune-datasets repo/repo --find-all --unstore --datasets icSrcSchema --no-confirm DECam/processing/*
echo Blasting logs
butler prune-datasets repo/repo --find-all --unstore --datasets *_log --no-confirm DECam/processing/*
echo Blasting deepCoadd_directWarp
butler prune-datasets repo/repo --find-all --unstore --datasets deepCoadd_directWarp --no-confirm DECam/processing/*
echo Blasting deepCoadd_psfMatchedWarp
butler prune-datasets repo/repo --find-all --unstore --datasets deepCoadd_psfMatchedWarp --no-confirm DECam/processing/*
echo Blasting directWarp
butler prune-datasets repo/repo --find-all --unstore --datasets directWarp* --no-confirm DECam/processing/*
echo Blasting overscanRaw
butler prune-datasets repo/repo --find-all --unstore --datasets overscanRaw --no-confirm DECam/processing/*
echo "Blasting coadds from metadetect (they're easy to recalculate)"
butler prune-datasets repo/repo --find-all --unstore --datasets deep_1p_Coadd --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_1m_Coadd --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_2p_Coadd --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_2m_Coadd --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_noshear_Coadd --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_1p_Coadd_calexp --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_1m_Coadd_calexp --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_2p_Coadd_calexp --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_2m_Coadd_calexp --no-confirm DECam/processing/*
butler prune-datasets repo/repo --find-all --unstore --datasets deep_noshear_Coadd_calexp --no-confirm DECam/processing/*

# the above steps leave the director structures behind, recursively delete to clear them (save the INODES)
rm -r repo/repo/DECam/processing/*/*/postISRCCD
rm -r repo/repo/DECam/processing/*/*/icSrc
rm -r repo/repo/DECam/processing/*/*/icExp
rm -r repo/repo/DECam/processing/*/*/icExpBackground
rm -r repo/repo/DECam/processing/*/*/icSrcSchema
rm -r repo/repo/DECam/processing/*/*/calexp
rm -r repo/repo/DECam/processing/*/*/*_log
rm -r repo/repo/DECam/processing/*/*/deepCoadd_directWarp
rm -r repo/repo/DECam/processing/*/*/deepCoadd_psfMatchedWarp
rm -r repo/repo/DECam/processing/*/*/directWarp*
rm -r repo/repo/DECam/processing/*/*/overscanRaw*
rm -r repo/repo/DECam/processing/meta_4a/*/deep_1p_Coadd*
rm -r repo/repo/DECam/processing/meta_4b/*/deep_1p_Coadd_calexp
rm -r repo/repo/DECam/processing/meta_4a/*/deep_1m_Coadd*
rm -r repo/repo/DECam/processing/meta_4b_1m/*/deep_1m_Coadd_calexp
rm -r repo/repo/DECam/processing/meta_4a/*/deep_2p_Coadd*
rm -r repo/repo/DECam/processing/meta_4b_2p/*/deep_2p_Coadd_calexp
rm -r repo/repo/DECam/processing/meta_4a/*/deep_2m_Coadd*
rm -r repo/repo/DECam/processing/meta_4b_2m/*/deep_2m_Coadd_calexp
rm -r repo/repo/DECam/processing/meta_4a/*/deep_noshear_Coadd*
rm -r repo/repo/DECam/processing/meta_4b_noshear/*/deep_noshear_Coadd_calexp

# zip-up the submit directory to save the INODES
zip -rm submit.zip submit

# this leaves the essential datasets in CLN/repo/repo/DECam/processing/calexp_* which should be archived (and used later in the future!)

