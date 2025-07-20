#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=200GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH -J select_visit_cluster_name
#SBATCH -o slurm_outputs/select_visit_cluster_name-%j.out
#SBATCH -e slurm_outputs/select_visit_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PY_SCRIPTS="py_scripts"
BAND=process_band
FWHM=fwhm_cut
ELLIP=ellip_cut
FWHM_r=fwhm_cut_r
ELLIP_r=ellip_cut_r

# navigate to .../A85
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# add the python_scripts from lovoccs_pipe to the PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:${PY_SCRIPTS}"

# by default this acts on all bands
# unfortunately this has to be run in serial since I cannot write to the repository in parallel

python python_scripts/check_visit/select_visit.py "${CLN}" "r" "${FWHM_r}" "${ELLIP_r}"

# y-band usually fails here, so skip it for now
for BAND in u g i z; do
	python python_scripts/check_visit/select_visit.py "${CLN}" "${BAND}" "${FWHM}" "${ELLIP}"
done


echo "Repairing headers..."

for BAND in u g i r z; do
	python python_scripts/misc/repair_calexp_observer.py "${CLN}" "repo/repo" "calexp" "DECam/processing/calexp_${BAND}"
done


echo "Creating the skymap..."

# create the skymap
#butler make-discrete-skymap -C configs/skymap_config.py --collections DECam/processing/quality_detectors_u,DECam/processing/quality_detectors_g,DECam/processing/quality_detectors_i,DECam/processing/quality_detectors_r,DECam/processing/quality_detectors_z --skymap-id "${CLN}_skymap" repo/repo DECam
butler register-skymap -C configs/skymap_config.py repo/repo


echo "Plotting the skymap..."

# plot the skymap
python python_scripts/render_data/draw_skymap_new.py "repo/repo" "${CLN}"

