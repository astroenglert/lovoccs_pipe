import sys
import os

from pathlib import Path

import numpy as np

from scipy.optimize import curve_fit
from scipy.integrate import cumulative_trapezoid as cumtrapz

import matplotlib.pyplot as pl
import matplotlib.ticker as ticker

from astropy.io import ascii
from astropy.table import Table, hstack, vstack
from astropy.cosmology import FlatLambdaCDM
from astropy.convolution import convolve, Gaussian1DKernel

# homebrew-imports here
# our default list of templates
from .bayesian_photo_z import large_sample_quality_check_plots, generic_quality_check_plots

from ..configs.photo_z_config import filter_map, error_tag, external_cache, mod_chi2_cut, odds_cut, prior_choice, truncate_bluest, sn_cuts

table_filename = "ALL-A2029_dered_dezp_gals_matched_specz_zphot.csv"
table = ascii.read(table_filename)

spec_z = table["z_spec"]
photo_z = table["z_phot"]

outliers_mask = photo_z > spec_z + 0.15 * (1.0 + spec_z)
outlier_table = table[outliers_mask]

large_sample_quality_check_plots(outlier_table,outlier_table,filter_map,z_spec_header="z_spec",output_filepath="./OUTLIER-A2029-")

# then delete zphot column
zphot_columns = [
    'z_phot', 'z_median', 'z_mode', 'z_ml', 
    'odds', 'mod_chi2', 
    'u_flux_at_z_phot', 'g_flux_at_z_phot', 'r_flux_at_z_phot', 'i_flux_at_z_phot', 'z_flux_at_z_phot'
]
outlier_table.remove_columns(zphot_columns)
outlier_table.write("OUTLIER-A2029_zphot_outliers.csv", format="ascii.csv", overwrite=True)



