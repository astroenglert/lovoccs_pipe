# lovoccs_pipe

A pipeline for processing DECam observations for the Local Volume Complete Cluster Survey (LoVoCCS).

This pipeline is built to walk the user through specific "processing steps" to go from raw images of galaxy clusters stored on the [NOIRLab Science Archive](https://astroarchive.noirlab.edu/) to lensing science. Going from raw data to coadds is carried out via a user's installation of the [LSST Science Pipelines](https://pipelines.lsst.io/) (LSP) and going from coadds/catalogs to lensing products is carried out through our own 'homebrewed' scripts. This build relies on LSP `v26.0.0` and relies on a modified version of the DECam Data Release Pipeline: `DRP-LoVoCCS.yaml`.

After making a clone of the repository somewhere, create a new folder where you want to process a cluster. Make a copy of the script `run_steps_Gen3.sh` and place it in this new folder, then update the variables at the top of the script as specified in the comments.

Running the pipeline is handled entirely by `run_steps_Gen3.sh`, which has a list of the 29 processing steps (and some notes) at the end of the file. After filling in the directories, resource allocation, and cluster name at the top of the script, you're ready to begin! To run a given processing step, comment out the line containing the command and run `bash run_steps_Gen3.sh` from the command line. `run_steps_Gen3.sh` uses text-replacement to fill lines in a template (see `processing_step_templates/`); the template is then saved to a folder specified in run_steps and automatically submitted to Slurm by calling `sbatch {PROCESSING_STEP}.sh`. The copy that is saved to a folder is ready to be manually submitted by the user via `sbatch` and is saved in case there are some manual tweaks that need to be made.


# Installation 

First, install LSP v26 somewhere (`~/lsst_stack_v26_0_0/`) and create a directory containing calibrations and refcats (`~/calib_catalog_repo`). Then download lovoccs_pipe in some directory (`~/processing`) via

``` git clone https://github.com/astroenglert/lovoccs_pipe.git ```

Create a new directory where you want to run the processing (e.g. A85)

``` mkdir ~/processing/A85/ ```

Copy `run_steps_Gen3.sh` into this new directory

``` cp ~/processing/lovoccs_pipe/run_steps_Gen3.sh ~/processing/A85/run_steps_Gen3.sh ```

Update the variables at the top of `run_steps_Gen3.sh` to reflect the cluster you want to process and the relevant directories

```
CLUSTER_NAME='A85'
...
NODES=10 
CORES=180 # multiple of NODES only!
RAM=1300 # multiple of NODES only!
WALL_TIME=48:00:00 
...
AUTO_PIPELINE_DIR="~/processing/lovoccs_pipe" # path to your installation of lovoccs_pipe
TEMPLATE_DIR="${AUTO_PIPELINE_DIR}/processing_step_templates" # shouldn't need to be adjusted
LOAD_PIPELINE_PATH="~/lsst_stack_v26_0_0/loadLSST.bash" # path to your v26.0.0 LSP installation
CLUSTER_DIR="~/processing/${CLUSTER_NAME}" # the directory in which you want to process
PROCESSING_STEP_DIR="${CLUSTER_DIR}/processing_step" # shouldn't need to be adjusted
CALIB_CATALOG_REPO="~/calib_catalog_repo" # location of calibrations and refcats
```

By default, the above directories and the filepaths below are all set for usage on Oscar/CCV here at Brown, but you can customize them as needed.

There are a handful of additional auxilliary files that we use throughout processing, these include:

1. A database containing calibrations and reference catalogs; the CALIB_CATALOG_REPO variable stores the path
   - By default, `ingest_data` expects to find a [reference catalog in the Gen3 format](https://pipelines.lsst.io/v/v26_0_0/modules/lsst.meas.algorithms/creating-a-reference-catalog.html) of type `refcat` located in `{CALIB_CATALOG_REPO}/gen3_formatted_new/{CLUSTER_NAME}/refcat`
   - Later steps also rely on the same set of catalogs stored as .csv's in `{CALIB_CATALOG_REPO}/catalogs_new/{CLUSTER_NAME}`
   - The calibrations, including `skyframe`, `bfk_kernel`, `flat`, `bias`, and `fringe` frames, should be [created w. the LSP Gen3](https://pipelines.lsst.io/v/v26_0_0/modules/lsst.cp.pipe/constructing-calibrations.html) and certified into a collection called `DECam/calib/certified`. This collection should be exported and saved in `{CALIB_CATALOG_REPO}/import_ready_calibs_6/`
   - These, and a list of all available refcats, also have filepaths specified in `python_scripts/configs/processing_step_configs.sh`

2. A directory containing colorterms and reference stars to be used during photometric calibration
   - The filepaths can be configured via the config-file `python_scripts/configs/photometric_correction_config.py`
   - Columns of the .csv are keyed by band and the rows specify the N-th order term in the polynomial (e.g. row-1 is the linear coefficient)
   - By default we calibrate using the g-i color of the corresponding reference catalog
   - The reference stars used to generate the terms are optional, but very useful to plot

3. A handful of other misc. files stored in `python_scripts/configs/processing_step_configs.sh`
   - SPECZ_DB is a database containing files used for reference by BPZ (just during plotting as a quality-check)
      - The expected format is a .csv with at least three headers: ra/dec/z
      - Should be titled by the cluster name, e.g. `{SPECZ_DB}/{CLUSTER_NAME}_ned_select.csv`; this can be switched by updating the photo_z processing step
   - EXT_DB contains extinction maps; by default we use IRAS maps
      - We use the default .fits file that can be (queried directly)[https://irsa.ipac.caltech.edu/applications/DUST/]
      - The expected directory structure is again based on the cluster name: `{EXT_DB}/{CLUSTER_NAME}.fits`
   - XRAY_DB specifies some archival X-Ray data that can be drawn over mass_maps for reference
      - Again, these are just .fits files saved in: `{XRAY_DB}/{CLUSTER_NAME}/Chandra/broad_flux_smoothed.fits`


# Usage

You're ready to process! Go to the bottom of the script, and uncomment the first step:

```
...
# == RUNNING STEPS == #

# == STEP0 NOTES == # 
# After copying and pasting run_steps_Gen3 into a directory with the Cluster Name, create_output creates a series of folders and python scripts. This only takes a few seconds to run

# == STEP0 COMMAND == #
create_output
...
```

And call run_steps from the command-line via `bash run_steps_Gen3.sh`. Each step is executed in exactly the same way, just don't forget to re-comment the previous step after verifying the step completed successfully!


# Citations
Works that make use of this pipeline should cite [Fu et al. 2021](https://doi.org/10.3847/1538-4357/ac68e8), [Englert et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025AAS...24541207E/abstract), and [Englert, Dell'Antonio, Montes (2025)](https://ui.adsabs.harvard.edu/abs/2025arXiv250523551E/abstract). Citations for LSP should also be included (see [Citing and acknowledging the LSST Science Pipelines](https://pipelines.lsst.io/#citing-and-acknowledging-the-lsst-science-pipelines)).


# Contributing

Please open an issue if there are major changes that need to be made, for small tweaks/bugfixes feel free to open a pull request.
