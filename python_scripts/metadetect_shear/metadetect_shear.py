
from typing import Any, ClassVar

import numpy as np

from numpy.fft import rfft2, irfft2, rfftfreq

import scipy as sp

import astropy.units as u

from metadetect.masking import (
    _build_square_apodization_mask,
    make_foreground_apodization_mask,
    make_foreground_bmask,
)

from metadetect.lsst.metacal_exposures import get_metacal_exps_fixnoise, get_metacal_exps

from lsst.daf.butler import Butler

import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
from lsst.afw.table import SimpleCatalog

from lsst.cp.pipe.utils import CovFastFourierTransform

from lsst.meas.algorithms import LoadReferenceObjectsConfig, ReferenceObjectLoader, getRefFluxField

from lsst.pex.config import ConfigField, Field, ChoiceField, ListField
from lsst.pipe.base import (
    InputQuantizedConnection,
    NoWorkFound,
    OutputQuantizedConnection,
    QuantumContext,
    PipelineTask,
    PipelineTaskConfig,
    PipelineTaskConnections,
    Struct,
)

import lsst.pipe.base.connectionTypes as cT

# I'll go with the pipetask approach for this, first I need to make the connections
# this will operate per-patch, per-band, hence the dimensions
class ShearCoaddConnections(PipelineTaskConnections,
                            dimensions=("tract","skymap","patch","band"),
                            defaultTemplates={"inputCoaddName": "deep"}
                           ):
    
    input_coadd = cT.Input(
        doc="Input coadd to mask and shear using metadetect",
        name="{inputCoaddName}Coadd",
        storageClass="ExposureF",
        dimensions=("tract","patch","skymap","band"),
    )
    
    ref_cat = cT.PrerequisiteInput(
        doc="Reference catalog used to mask bright objects",
        name="gaia",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )
    
    coadd_noshear = cT.Output(
        doc="noshear coadd",
        name="{inputCoaddName}_noshear_Coadd",
        storageClass="ExposureF",
        dimensions=("tract","patch","skymap","band"),
    )
    
    coadd_1p = cT.Output(
        doc="1p sheared coadd",
        name="{inputCoaddName}_1p_Coadd",
        storageClass="ExposureF",
        dimensions=("tract","patch","skymap","band"),
    )
    
    coadd_1m = cT.Output(
        doc="1m sheared coadd",
        name="{inputCoaddName}_1m_Coadd",
        storageClass="ExposureF",
        dimensions=("tract","patch","skymap","band"),
    )
    
    coadd_2p = cT.Output(
        doc="2p sheared coadd",
        name="{inputCoaddName}_2p_Coadd",
        storageClass="ExposureF",
        dimensions=("tract","patch","skymap","band"),
    )
    
    coadd_2m = cT.Output(
        doc="2m sheared coadd",
        name="{inputCoaddName}_2m_Coadd",
        storageClass="ExposureF",
        dimensions=("tract","patch","skymap","band"),
    )
    
    def __init__(self, *, config=None):
        super().__init__(config=config)
        
        # if you don't want to mask bright stars, remove the ref_cat from prereq
        if not self.config.doMaskBrightStars:
            self.prerequisiteInputs.remove("ref_cat")
    

class ShearCoaddConfig(PipelineTaskConfig, pipelineConnections=ShearCoaddConnections):
    
    """ Configuration definition for ShearCoaddTask """
    
    doMaskBrightStars = Field(
        dtype=bool,
        default=True,
        doc="Mask bright stars before shearing?",
    )
    
    ref_loader = ConfigField(
        dtype=LoadReferenceObjectsConfig,
        doc="Reference object loader used for bright-object masking.",
    )
    
    noise_image = ChoiceField(
        dtype=str,
        doc="Method used to generate noise image, options are 'gaussian_field' or 'bootstrap' ",
        default="gaussian_field",
        allowed={
            "gaussian_field" : "generate a gaussian field with the same spectrum",
            "bootstrap" : "generate a field by resampling bkg-pixels with replacement",
        }
    )
    
    grow_detections = Field(
        dtype=int,
        doc="Pixels by which to grow detection masks before measuring the noise",
        default=16,
    )
    
    seed_rng = Field(
        dtype=bool,
        doc="Seed the rng with a specific value?",
        default=False,
    )
    
    seed_rng_val = Field(
        dtype=int,
        doc="Seed used to initiate the rng, ignored if seed_rng disabled",
        default=1,
    )
    
    bright_star_filter = Field(
        dtype=str,
        doc="Flux field in ref_cat used for building the bright star mask",
        default="phot_g_mean_mag",
    )
    
    do_gaussian_field = Field(
        dtype=bool,
        doc="Create a gaussian field w. equivalent power spectrum? Otherwise a simple noise is generated. Ignored of do_noise is False",
        default=True,
    )
    
    do_noise = Field(
        dtype=bool,
        doc="Create a noise image?",
        default=True,
    )
    


# I'm using the metadetection_shear from tickets/DM-40513 as a base for this
class ShearCoaddTask(PipelineTask):
    """ A PipelineTask that applies four shears to coadds following the algorithm outlined in Sheldon+23 """
    
    _DefaultName: ClassVar[str] = "shearCoadd"
    ConfigClass: ClassVar[type[ShearCoaddConfig]] = ShearCoaddConfig
    
    config = ShearCoaddConfig()
    
    # is this necessary? Only if I need to tweak something on initializing the object
    def __init__(self, *, initInputs: dict[str, Any] | None = None, **kwargs: Any):
        super().__init__(initInputs=initInputs, **kwargs)
        
    # same with runQuantum, I don't need to tweak it unless there is something essential I need to override
    def runQuantum(self,
                   butlerQC: QuantumContext,
                   inputRefs: InputQuantizedConnection,
                   outputRefs: OutputQuantizedConnection,
                  ):
        
        # load the reference catalog
        ref_loader = ReferenceObjectLoader(
            dataIds=[ref.datasetRef.dataId for ref in inputRefs.ref_cat],
            refCats=[butlerQC.get(ref) for ref in inputRefs.ref_cat],
            name=self.config.connections.ref_cat,
            config=self.config.ref_loader,
            log=self.log,
        )
        
        ref_cat = ref_loader.loadRegion(butlerQC.quantum.dataId.region,filterName=self.config.bright_star_filter)
        
        # load the exposure
        shear_me = butlerQC.get(inputRefs.input_coadd)
        
        # seed the rng
        if not self.config.seed_rng:
            # pull numbers from the cluster_name, stored in the skymap
            # then append the patch number to that sequence
            # and append a number for each band
            # simple and reproducible, this doesn't need to be anything more than that!
            coadd_id = inputRefs.input_coadd.dataId
            cluster_numbers = ''.join(num for num in coadd_id['skymap'] if num.isdigit())
            patch_number = str( coadd_id['patch'] )
            band_number = str( "ugriz".index(shear_me.filter.bandLabel) )
            seed = int(cluster_numbers + patch_number + band_number)
            gen = np.random.default_rng(seed = seed)
        else:
            gen = np.random.default_rng(seed = self.config.seed_rng_val)
        
        # now that everything is loaded, start actually running the task!
        outputs = self.run(coadd=shear_me,ref_cat=ref_cat.refCat,generator=gen)
        
        butlerQC.put(outputs, outputRefs)
    
    
    def run(
        self,
        coadd: afwImage.ExposureF,
        ref_cat: SimpleCatalog,
        generator: np.random.Generator
    ) -> Struct:
        """
        Using the algorithm from metadetect, shear each patch.
        
        Parameters
        ----------
        coadd : `~lsst.afw.image.ExposureF`
            A coadd to shear
        generator: np.random.Generator
            A generator to use for creating the noise image
        ref_Cat : `lsst.afw.table.SimpleCatalog`
            A reference catalog to use while masking bright stars
        
        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Structure conaining the output sheared coadds:
            
            ``coadd_noshear``
                coadd with no shear applied
            ``coadd_1p``
                coadd with +0.01 shear along e1
            ``coadd_2p``
                coadd with +0.01 shear along e2
            ``coadd_1m``
                coadd with -0.01 shear along e1
            ``coadd_2m``
                coadd with -0.01 shear along e2
                
        """
        
        # first grow the detection mask (growMasksIsotropic)
        mask = coadd.getMask().clone()
        growMasksIsotropic(mask, radius=self.config.grow_detections, maskNameList=['DETECTED'], maskValue="DETECTED")
        
        # second create the noise image (covariance_fourier and realize_field, add alternate for bootstrap)
        # load the detection mask
        full_mask = mask
        mask_array = full_mask.array >> int(np.log2(full_mask.getPlaneBitMask('DETECTED'))) & 1
        
        if self.config.do_noise & self.config.do_gaussian_field:
            
            # now compute the covariance
            cov, cov_err = covariance_fourier(coadd.image.array, 1 - mask_array, maxRangeCov=300)
            noise_array = realize_field(cov, size=4200, gen=generator)
            
            # scale the simulated noise to match the variance-plane?
            noise_array = noise_array * np.median(coadd.variance.array)/np.std(noise_array)
            
            # create a copy of the coadd image to store the noise_image; copy over metadata
            noise_image = coadd.clone()
            noise_image.image.array[:,:] = noise_array
        
        elif self.config.do_noise:
            
            # if gaussian field is disabled, do a simple model
            noise_image = simple_noise(coadd, 1 - mask_array, gen=generator)
        
        else:
            
            # if None, create_shear_exp catches this and falls back to no noise_image
            noise_image = None
        
        # third apodize-mask edges (apply_apodized_edge_mask)
        apply_apodized_edge_mask(coadd, noise_image, AP_RAD=1.5)
        
        # fourth select and apodize-mask bright stars (apply_apodized_bright_mask and collect_bright_stars)
        bright_info = collect_bright_stars(ref_cat, filter_name='phot_g_mean_mag', min_radius=5)
        apply_apodized_bright_mask(coadd, noise_image, bright_info, EXPAND_RAD=16, AP_RAD=1.5)
        
        # one last thing, clear the detection masks since detection will be re-run!
        # shamelessly taken from lsst.meas.algorithms.detection
        # not needed, detection clears det mask before re-detecting
        #mask = coadd.getMask()
        #mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))
        
        # fifth run the shear algorithm from metadetect (shear_exp)
        shear_exp, shear_noise_exp = create_shear_exp(coadd, noise_image)

        # add original PSF back to metadetect coadds
        final_exp = {}
        for s in shear_exp.keys():
            c = afwImage.ExposureF(src=shear_exp[s],deep=True)
            c.setPsf(coadd.getPsf())
            final_exp[s] = c
        
        # lastly return the outputs
        return Struct(
            coadd_noshear = final_exp['noshear'],
            coadd_1p = final_exp['1p'],
            coadd_2p = final_exp['2p'],
            coadd_1m = final_exp['1m'],
            coadd_2m = final_exp['2m'],
        )
        
        
    
# function to actually apply a shear and collect the four new exposures
def create_shear_exp(exp, noise_exp, types=('noshear', '1p', '1m', '2p', '2m')):
    """
    Function to shear the exposure using tools from metadetect.lsst
    
    Parameters
    ----------
    exp: ExposureF
        Exposure to apply a shear to
    noise_exp: ExposureF
        Noise image with the same spectral characterstics as exp
    types: array
        Array containing the types of shears to apply
    
    Returns
    -------
    shear_dict: dict
        Dictionary containing the five image (shear/no-shear)
    shear_noise_dict: dict
        Dictionary containing the five sheared noise images
    
    """
    
    if noise_exp == None:
        shear_dict = get_metacal_exps(exp, types=types)
        shear_noise_dict = None
    else:
        shear_dict, shear_noise_dict = get_metacal_exps_fixnoise(exp, noise_exp, types=types)
    
    
    return shear_dict, shear_noise_dict
    

def simple_noise(exp, mask_array, gen):
    """
    Helper function to run a simple noise model
    
    Parameters
    ----------
    exp: Exposure
        Exposure to draw noise from
    mask_array: numpy.array
        Mask, 1/True are background pixels
    gen: Generator
        RNG
    
    Returns
    -------
    noise_exp: Exposure
        Noise exposure with metadata inherited from exp
        
    """
    
    noise_exp = exp.clone()
    
    # this is correct
    noise_exp.image.array[:,:] = gen.normal(scale=np.sqrt(np.median(exp.variance.array)), size=noise_exp.image.array.shape)
    
    return noise_exp
    

# there are more masks to consider, e.g. nearby galaxies, globlars, dwarfs, etc.
# but lovoccs has a carefully curated set of fields, so this shouldn't be a huge issue and can be handled on a case-by-case basis other than the base star-mask
def collect_bright_stars(ref_cat, filter_name='phot_g_mean_mag', min_radius=5):
    """
    Collect bright stars to mask from a ref_cat
    
    Parameters
    ----------
    ref_cat: SimpleCatalog
        reference catalog containing stars in the field
    filter_name: string
        name of the filter in ref_cat to use
    min_radius: float
        minimum radius to mask around each star (in arcseconds)
        
    Returns
    -------
    bright_info: dict
        Dictionary containing arrays storing the 'ra', 'dec' of bright stars and the 'radius_pixels' to mask surrounding them
    
    """
    
    # Radius from DES Y6 results, we may need to be more agressive? Our data is deeper so the stars will be a touch fatter
    radius = lambda g: 10**( (0.004432 * g**2) - (0.2257 * g) + 2.996 ) # in arcseconds
    
    # convert nJy to data; collect the flux field using the table schema and convert it to AB-mag
    flux_field = getRefFluxField(ref_cat.schema,filter_name)
    star_mag = u.Quantity(ref_cat[flux_field], u.nJy).to_value(u.ABmag)
    
    # assign the radius to each star
    # need to load the pixel-size properly, later problem!
    radius_arcsec = radius(star_mag)
    radius_arcsec[radius_arcsec <= 5] = 5
    radius_pixels = radius_arcsec/0.263
    
    # collect the ra/dec and convert them to degrees
    coord_key = ref_cat.getCoordKey()
    ra = u.Quantity(ref_cat[coord_key.getRa()],u.radian).to_value(u.degree)
    dec = u.Quantity(ref_cat[coord_key.getDec()],u.radian).to_value(u.degree)
    bright_info = dict(ra=ra, dec=dec, radius_pixels=radius_pixels)
    
    return bright_info


# in hindsight... I should've just wrapped the single exp/noise_exp into a mbexp object rather than re-writing all these functions. Oh well, it was a good learning experience at least :D
# modifed version of apply_apodized_bright_masks_mbexp built to run on a single exposure
def apply_apodized_bright_mask(exp, noise_exp, bright_info, EXPAND_RAD=16, AP_RAD=1.5):
    """
    Apply an apodized mask to bright stars
    
    Parameters
    ----------
    exp: ExposureF
        Exposure to apply an apodized mask to
    noise_exp: ExposureF
        Noise exposure, which is given the same mask
    bright_info: dict
        Dictionary containing 'ra', 'dec', and 'pixel_radius' to mask
    EXPAND_RAD: float
        Radius within which all detections are ignored
    AP_RAD: float
        Parameter for the apodization kernel
    
    """
    
    # create new mask planes for BRIGHT and BRIGHT_EXPANDED masks
    afwImage.Mask.addMaskPlane('BRIGHT')
    afwImage.Mask.addMaskPlane('BRIGHT_EXPANDED')
    bright = afwImage.Mask.getPlaneBitMask('BRIGHT')
    bright_expanded = afwImage.Mask.getPlaneBitMask('BRIGHT_EXPANDED')
    bad = afwImage.Mask.getPlaneBitMask('BAD')
    
    # collect the wcs and origin of the exp coordinate system
    wcs = exp.getWcs()
    xy0 = exp.getXY0()
    
    # collect the coordinates of bright objects to mask and other prereqs for the masks
    xm, ym = wcs.skyToPixelArray(ra=bright_info['ra'],dec=bright_info['dec'], degrees=True)
    xm -= xy0.x
    ym -= xy0.y
    rm = bright_info['radius_pixels']
    dims = exp.image.array.shape
    
    # now actually build the mask
    ap_mask = make_foreground_apodization_mask(
        xm=xm,
        ym=ym,
        rm=rm,
        dims=dims,
        symmetrize=False,
        ap_rad=AP_RAD,
    )
    
    # apply to the image/mask
    msk = np.where(ap_mask < 1)
    
    exp.image.array[msk] *= ap_mask[msk]
    exp.variance.array[msk] = 1e8 #np.inf
    exp.mask.array[msk] |= bright
    exp.mask.array[msk] |= bad
    
    if noise_exp != None:
        noise_exp.image.array[msk] *= ap_mask[msk]
        noise_exp.variance.array[msk] = 1e8 #np.inf
        noise_exp.mask.array[msk] |= bright
        noise_exp.mask.array[msk] |= bad
    
    # build the expanded mask as well
    expanded_bmask = make_foreground_bmask(
        xm=xm,
        ym=ym,
        rm=rm + EXPAND_RAD,
        dims=dims,
        symmetrize=False,
        mask_bit_val=bright_expanded,
    )
    
    # this is computationally cheap, lazy solution is to just change the mask_bit to bad 
    expanded_bmask_bad = make_foreground_bmask(
        xm=xm,
        ym=ym,
        rm=rm + EXPAND_RAD,
        dims=dims,
        symmetrize=False,
        mask_bit_val=bad,
    )
    
    # apply the extended mask
    exp.mask.array[:,:] |= expanded_bmask
    exp.mask.array[:,:] |= expanded_bmask_bad
    if noise_exp != None:
        noise_exp.mask.array[:,:] |= expanded_bmask
        noise_exp.mask.array[:,:] |= expanded_bmask_bad
    
    return


# modifed version of apply_apodized_edge_masks_mbexp built to run on a single exposure
def apply_apodized_edge_mask(exp, noise_exp, AP_RAD=1.5):
    """
    Mask the edges of an exp and its corresponding noise image
    
    Parameters
    ----------
    exp: lsst.afw.image.Exposure
        The data to mask. The image and mask are modified.
    noise_exp: lsst.afw.image.Exposure
        The noise image to mask. Again, the image and mask are modified
    AP_RAD: float
        Parameter passed to the apodiz kernel
    
    """
    # create new mask planes to ignore during detection later
    afwImage.Mask.addMaskPlane('APODIZED_EDGE')
    edge = afwImage.Mask.getPlaneBitMask('APODIZED_EDGE')
    bad = afwImage.Mask.getPlaneBitMask('BAD')
    
    # now create the actual mask
    ap_mask = np.ones_like(exp.image.array)
    _build_square_apodization_mask(AP_RAD, ap_mask)
    
    # update the variance/image/masks of exp and noise_exp
    msk = np.where(ap_mask < 1)
    exp.image.array[msk] *= ap_mask[msk]
    exp.variance.array[msk] = 1e8 #np.inf
    exp.mask.array[msk] |= edge
    exp.mask.array[msk] |= bad
    
    if noise_exp != None:
        noise_exp.image.array[msk] *= ap_mask[msk]
        noise_exp.variance.array[msk] = 1e8 #np.inf
        noise_exp.mask.array[msk] |= edge
        noise_exp.mask.array[msk] |= bad
    
    return


# a modifed version of growMasks from lsst.ip.isr which grows the masks isotropically rather than manhattan
def growMasksIsotropic(mask, radius=16, maskNameList=['DETECTED'], maskValue="DETECTED"):
    """Grow a mask isotropically by an amount and add to the requested plane.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        Mask image to process.
    radius : scalar
        Amount to grow the mask.
    maskNameList : `str` or `list` [`str`]
        Mask names that should be grown.
    maskValue : `str`
        Mask plane to assign the newly masked pixels to.
    """
    
    if radius > 0:
        thresh = afwDetection.Threshold(mask.getPlaneBitMask(maskNameList), afwDetection.Threshold.BITMASK)
        fpSet = afwDetection.FootprintSet(mask, thresh)
        fpSet = afwDetection.FootprintSet(fpSet, rGrow=radius, isotropic=True)
        fpSet.setMask(mask, maskValue)
    
    return


# function to run the fourier-space calculation of covariance
# uses method from Appendix A Astier+19 and borrows from the implementation of it in cp_pipe
# assumes quadrant-symmetry, good enough for the noise based on earlier tests
# propagate detection mask out 16px, 4xFWHM, to get the noise spectrum
# beyond 300px lags we start to lose the "covariance signal" (Cij / (V/root(N))) due to a lack of samples
def covariance_fourier(image, mask, maxRangeCov=300, fftShape=None):
    """
    Parameters:
      image: numpy array; a numpy array to compute the covariance for
      mask: numpy array; a numpy array containing pixels to mask (1 = masked)
      maxRangeCov: int; the maximum covariance calculated
      fftShape: numpy array; Optional: the size of the FFT, otherwise uses LSST defaults
    
    Returns:
      cov_array: numpy array; the covariance of an imge
      cov_err_array: numpy array; the error in covariance of an image
      
    """
    
    # if this isn't specified, follow the default from lsst.cp.pipe.ptc.cpPtcExtract
    if fftShape == None:
        im_shape = np.array(image.shape)
        s = im_shape + maxRangeCov
        tempSize = np.array(np.log(s)/np.log(2.)).astype(int)
        fftSize = np.array(2**(tempSize+1)).astype(int)
        fftShape = (fftSize[0],fftSize[1])
    
    c = CovFastFourierTransform(image, mask, fftShape, maxRangeCov)
    
    try:
        covDiffAstier = c.reportCovFastFourierTransform(maxRangeCov)
    except ValueError:
        print('Not enough pixels covering the requested covariance range in x/y')
        return np.nan, np.nan
    
    # now format the output nicely
    cov_array = np.zeros((maxRangeCov+1,maxRangeCov+1))
    cov_err_array = np.zeros((maxRangeCov+1,maxRangeCov+1))
    
    for entry in covDiffAstier:
        cov_array[entry[0],entry[1]] = entry[3]
        cov_err_array[entry[0],entry[1]] = entry[2]/np.sqrt(entry[4])
    
    return cov_array, cov_err_array


# from a covariance, draw a random field with the same power spectrum
# pads cov_array with zeros to match size+1, e.g. assumes beyond cov_array field is uncorrelated
def realize_field(cov_array, size=4000, gen=None):
    """
    Parameters:
      cov_array: numpy array; an array containing the i,j-th pixel covariance (quadrant sym)
      size: int; an int specifying the size of the field to draw; should be even
      gen: Generator; generator for rng
    
    Returns:
      real_field: numpy array; an array storing the noise field
    
    """
    
    if gen == None:
        gen = np.random.default_rng()
    
    # stitch together the cov_array to create the full covariance
    cov_full = np.block([[np.flipud(np.fliplr(cov_array[:,1:])),np.flipud(cov_array)],[np.fliplr((cov_array[1:,1:])),(cov_array[1:,:])]])
    cov_shape = int( (cov_full.shape[0] - 1)/2 )
    
    # zero-pad covariance to center it, noise is uncorrelated on large scales
    temp = np.zeros((size+1,size+1))
    temp[int(size/2 - cov_shape):int(size/2 + cov_shape + 1),int(size/2 - cov_shape):int(size/2 + cov_shape + 1)] = cov_full
    
    # compute the fourier transform
    pk_full = rfft2(temp)
    
    # assemble random phases to offset each mode
    random_phase = np.exp(1j * gen.uniform(high=2*np.pi,size=pk_full.shape) )
    fourier_field = np.sqrt(pk_full) * random_phase
    real_field = ( irfft2(fourier_field) * fourier_field.shape[0] )[:-1,:]
    
    return real_field










