#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -J check_visit_process_band_cluster_name
#SBATCH -o slurm_outputs/check_visit_process_band_cluster_name-%j.out
#SBATCH -e slurm_outputs/check_visit_process_band_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
BAND=process_band
FWHM=fwhm_cut
ELLIP=ellip_cut

# navigate to .../A85
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# prevent implicit multithreading (otherwise tasks compete for resources)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p check_visit/band_summary_plots
mkdir -p check_visit/plots
mkdir -p check_visit/summary_tables

python -u python_scripts/check_visit/check_visit.py "${BAND}" "20"

python -u python_scripts/check_visit/draw_psf.py "${BAND}" "20" "${FWHM}" "${ELLIP}"


