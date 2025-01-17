This repository contains open-source python software that provides the weight and calibration factors associated with the first-year (S16A) shape catalog from the HSC SSP survey ([catalog description](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/), [catalog paper](https://ui.adsabs.harvard.edu/abs/2018PASJ...70S..25M/abstract), [simulation paper](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.3170M/abstract)).  The software is available for use under a BSD 3-clause license contained in the file [LICENSE.md](LICENSE.md).

Currently the script `gen_hsc_calibrations.py` is set up to read in a FITS catalog containing the outputs from Rubin’s LSST Science Pipelines, including re-Gaussianization shape measurements and cmodel magnitudes and errors, and produce a catalog with the additional quantities used for ensemble weak lensing shear estimation.  (The contents of `utilities.py` and of `data/` are accessed by `gen_hsc_calibrations.py`, but users need not interact with them directly.)  For information about dependencies, command-line arguments, and memory usage, please read the docstring for `gen_hsc_calibrations.py`; the docstring will also be output to the terminal if you call the script without any arguments.  Note that the script could easily be modified to serve as a python function that takes in and returns datasets while executing a longer workflow.  For instructions on how to use the outputs of this script, please see appendix A3 of the [catalog paper](https://ui.adsabs.harvard.edu/abs/2018PASJ...70S..25M/abstract).

# Limitations

This script implements the calibrations described in the [simulation paper](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.3170M/abstract), and comes with the following caveats:

1. There are a number of caveats associated with shear calibration via the particular simulations used for that purpose, as outlined in sections 5.4-5.7 and 6 of that paper.

2. The calibrations and weight functions were derived for samples that used the cuts outlined in section 5.1 or table 4 of the [catalog paper](https://ui.adsabs.harvard.edu/abs/2018PASJ...70S..25M/abstract).  Application of the routines to samples that fail those cuts is not recommended and the results could have unexpected features.  Similarly, use of different weights means that the corrections for weight bias are no longer applicable, and new ones would have to be derived.

3. The calibrations and weight functions were derived for the HSC S16A shear catalog.  Application to datasets with very different seeing characteristics or depth may not be appropriate, due to differences in blending-related biases and other factors. Users are responsible for the use of calibrations derived for their own data set, including estimating systematic uncertainties from doing so. 

# Communicating with the developers

If you have questions about this software, please open an [issue in this repository](https://github.com/PrincetonUniversity/hsc-y1-shear-calib/issues) to communicate with the developers.

# Attribution

If you use this software for work resulting in a publication, please do the following to provide attribution:

* Cite the [catalog paper](https://ui.adsabs.harvard.edu/abs/2018PASJ...70S..25M/abstract) and [simulation paper](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.3170M/abstract), and
* Use the HSC SSP acknowledgments given [here](https://hsc-release.mtk.nao.ac.jp/doc/index.php/acknowledging-hsc-2/).
