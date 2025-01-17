#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH -J download_noao_cluster_name
#SBATCH -o slurm_outputs/download_noao_cluster_name-%j.out
#SBATCH -e slurm_outputs/download_noao_cluster_name-%j.err

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name" 
LOAD_LSST="load_pipeline_path" 
CLUSTER_DIR="cluster_dir" 

# navigate to .../A85
cd ${CLUSTER_DIR}
mkdir -p raws

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# clear previous download tasks (if they exist)
rm raws/download_task_*.sh

# submit the noao query via python, save the download commands to download_raws.sh
echo "Querying Noirlab..."
python -u "python_scripts/noirlab_download/noao_query_raw.py" ${CLN}

python -u "python_scripts/noirlab_download/download_manager.py"

rm raws/download_task_*.sh
