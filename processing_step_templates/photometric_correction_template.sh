#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -J photometric_correction_cluster_name
#SBATCH -o slurm_outputs/photometric_correction_cluster_name-%j.out 
#SBATCH -e slurm_outputs/photometric_correction_cluster_name-%j.out 

# defining variables
# TEXT REPLACED FOR TEMPLATES
CLN="cluster_name"
LOAD_LSST="load_pipeline_path"
CLUSTER_DIR="cluster_dir" # UPDATE LATER
PHOTOM="photom_ref" # the refcat used for photometric_correction, e.g. ps1/sdss/legacy/des
PHOTOM_TAG="photom_dr" # tag specifying information about the data-release of the photom-refcat
#TODO Pass this via run_steps rather than hard-coded
CALIB_CATALOG_REPO="/gpfs/data/idellant/Clusters/calib_catalog_repo" 

# navigate to .../A85
cd ${CLUSTER_DIR}

# initalize the LSP (LSST Science Pipeline)
source ${LOAD_LSST}
setup lsst_distrib

# create an output directory for photometric_correction
mkdir photometric_correction_output


# first, de-redden catalog and refcat
# TODO how the filepaths here should come from a config rather than being hard-coded here
EXT_IMAGE="/gpfs/data/idellant/Clusters/galactic_extinction_database/data_extinction_irsa/${CLN}.fits"
INPUT_CAT="read_catalog_all_output/${CLN}_00-1111_all.csv"
REF_CAT="${CALIB_CATALOG_REPO}/catalogs_new/${CLN}/${PHOTOM}_${PHOTOM_TAG}_${CLN}.csv"

echo "Running the extinction correction!"
python python_scripts/photometric_correction/extinction_correction.py ${EXT_IMAGE} ${INPUT_CAT} "photometric_correction_output/${CLN}_dered.csv" "decam"
python python_scripts/photometric_correction/extinction_correction.py ${EXT_IMAGE} ${REF_CAT} "photometric_correction_output/${CLN}_${PHOTOM}_dered.csv" "${PHOTOM}"


# second, separate the stars and galaxies in our observations
echo "Separating stars and galaxies!"
python python_scripts/photometric_correction/separate_stars_galaxies.py "photometric_correction_output/${CLN}_dered.csv" "photometric_correction_output/${CLN}_dered_stars.csv" "photometric_correction_output/${CLN}_dered_gals.csv"


# third, matched the de-reddned refcat with the de-reddned star-catalog
echo "Matching the refcat!"
python python_scripts/misc/match_catalogs.py "photometric_correction_output/${CLN}_${PHOTOM}_dered.csv" "${PHOTOM}" "_ref" "photometric_correction_output/${CLN}_dered_stars.csv" "decam" "_cat" "photometric_correction_output/${CLN}_dered_stars_matched_${PHOTOM}.csv" "0.2"


# fourth, apply color_terms and fit for a linear-bias
# implicitly, scripts assume _cat and _ref tags are applied during matching!
echo "Comparing the catalogs (color_terms)"
python python_scripts/photometric_correction/color_terms.py "photometric_correction_output/${CLN}_dered_stars_matched_${PHOTOM}.csv" "decam" "${PHOTOM}" "photometric_correction_output/${CLN}_matched_residuals.csv"


# fifth, run the stellar_locus to get a linear-bias for the u-band. Uses a de-reddened star-catalog and creates a new csv w. the zp for u-band; saved to photom_ref_matched_residuals_stellar_locus.csv
echo "Running the locus-correction"
python python_scripts/photometric_correction/stellar_locus.py "photometric_correction_output/${CLN}_dered_stars.csv" "decam" "photometric_correction_output/${CLN}_matched_residuals.csv"


# sixth, apply the zp-corrections to the entire catalog
echo "Applying the zp-correction"
python python_scripts/photometric_correction/zero_point.py "photometric_correction_output/${CLN}_dered_stars.csv" "photometric_correction_output/${CLN}_matched_residuals.csv" "photometric_correction_output/${CLN}_matched_residuals_stellar_locus.csv" "photometric_correction_output/${CLN}_dered_dezp_stars.csv"

python python_scripts/photometric_correction/zero_point.py "photometric_correction_output/${CLN}_dered_gals.csv" "photometric_correction_output/${CLN}_matched_residuals.csv" "photometric_correction_output/${CLN}_matched_residuals_stellar_locus.csv" "photometric_correction_output/${CLN}_dered_dezp_gals.csv"

python python_scripts/photometric_correction/zero_point.py "photometric_correction_output/${CLN}_dered.csv" "photometric_correction_output/${CLN}_matched_residuals.csv" "photometric_correction_output/${CLN}_matched_residuals_stellar_locus.csv" "photometric_correction_output/${CLN}_dered_dezp.csv"


# and lastly, render figures comparing our corrections with other refcats
# this step is a little more painful, since we need to search for the existing refcats
#TODO somehow this should really be wrapped in a config-option, e.g. some list of refcats should be stored
# but until we implement that, a for loop over the possible catalogs and their corresponding instruments should do
POSSIBLE_CATALOGS=("des_dr2" "legacy_survey_dr9" "ps1_dr1" "sm_dr1" "sm_dr2")
INSTRUMENTS=("des" "legacy" "ps1" "sm" "sm")

for DEX in ${!POSSIBLE_CATALOGS[@]}; do
    
    # load the catalog and corresponding instrument
    #TODO these will need to be implemented with the IO
    CAT=${POSSIBLE_CATALOGS[$DEX]}
    CAT_PATH="${CALIB_CATALOG_REPO}/catalogs_new/${CLN}/${CAT}_${CLN}.csv"
    INS=${INSTRUMENTS[$DEX]}
    
    echo "Now comparing with ${CAT}"
    # if the catalog doesn't exist, continue the loop
    #TODO something is wrong with the line below, not sure what?
    #if [ ! -f ${CAT_PATH} ]; then continue; fi
    
    # dered and match each refcat 
    python python_scripts/photometric_correction/extinction_correction.py ${EXT_IMAGE} ${CAT_PATH} "photometric_correction_output/${CAT}_dered.csv" "${INS}"
    
    python python_scripts/misc/match_catalogs.py "photometric_correction_output/${CAT}_dered.csv" "${INS}" "_ref" "photometric_correction_output/${CLN}_dered_dezp_stars.csv" "decam" "_cat" "photometric_correction_output/${CLN}_dered_stars_matched_${CAT}.csv" "0.2"
    
    # now run the comparison
    python python_scripts/photometric_correction/compare_catalogs.py "photometric_correction_output/${CLN}_dered_stars_matched_${CAT}.csv" "decam" "${INS}"

done
