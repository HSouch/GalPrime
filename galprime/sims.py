import galprime
import numpy as np
import os

from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel

from numpy.random import randint
from scipy.signal import convolve2d

import galprime

from pebble import ProcessPool
from concurrent.futures import TimeoutError

from tqdm import tqdm

class GalPrimeError(Exception):
    pass


class GPrime():
    """ Class to create and run GalPrime simulations """
    
    def __init__(self, config=None, psfs=None, backgrounds=None, mag_kde=None):

        if config is None:
            self.config = galprime.default_config_params()
        elif type(config) == str:
            self.config = galprime.load_config_file(config)
        elif type(config) == dict:
            self.config = config
        else:
            raise GalPrimeError("Uncertain format for config file. File needs to be dict object or str filename.")
            sys.exit(1)
        
        self.binned_objects = None
        self.psfs = psfs
        self.backgrounds = backgrounds
        self.mag_kde = mag_kde

    def load_bins(self):
        self.binned_objects = galprime.bin_catalog(self.config)
    
    
    def pipeline(self, max_bins=None, mag_kde=None, process_method=None,
                progress_bar=False, debug=False, 
                 table_names=["ID","IMG", "X", "Y", "RA", "DEC", "MAGS", "R50S", "NS", "ELLIPS", "MASS_MED",
                              "ZPHOT", "sfProb", "I_R50", "PA", "R50_PIX", "BG_MEAN", "BG_MED", "BG_STD",
                              "NSEG_BGA", "CENT_BGA", "N_MASKED_BGA", "P_MASKED_BGA", "NSEG_BGSUB", "CENT_BGSUB", 
                              "N_MASKED_BGSUB", "P_MASKED_BGSUB"]):
        
        # Bin out catalog based on the catalog in the configuration file
        binned_objects = galprime.bin_catalog(self.config)
        max_bins = len(binned_objects.bins) if max_bins is None else max_bins
        verbose = self.config["VERBOSE"]

        if verbose:
            print("Running GalPRIME, version", galprime.__version__)
            print(len(binned_objects.bins), "in total. Running GalPRIME on", max_bins)


        # Create all of our necessary filestructure
        if not os.path.isdir(self.config["OUT_DIR"]):
            os.mkdir(self.config["OUT_DIR"])

        # Generate the full filestructure, including containers for the individual profiles, and the medians
        galprime.gen_filestructure(self.config["OUT_DIR"])
        bare_dir = self.config["OUT_DIR"] + "bare_profiles/"
        bgadded_dir = self.config["OUT_DIR"] + "bgadded_profiles/"
        bgsub_dir = self.config["OUT_DIR"] + "bgsub_profiles/"
        
        bare_median_dir = self.config["OUT_DIR"] + "bare_medians/"
        bgadded_median_dir = self.config["OUT_DIR"] + "bgadded_medians/"
        bgsub_median_dir = self.config["OUT_DIR"] + "bgsub_medians/"

        additional_dir = self.config["OUT_DIR"] + "additional_info/"

        
        # Run through bins and process using the method
        for i in range(max_bins):
            current_bin = binned_objects.bins[i]
            
            columns = current_bin.columns()
            
            # Generate the object KDE that will be used by this bin
            kde = galprime.object_kde(columns)
                        
            job_list, job_results = [], []
            # Set up the job pool
            print("Generating", self.config["N_MODELS"], "models for bin:", current_bin.bin_params)
            with ProcessPool(max_workers=self.config["CORES"]) as pool:
                for i in range(self.config["N_MODELS"]):
                    job_list.append(pool.schedule(process_method,
                                                args=(self, current_bin, kde),
                                                timeout=self.config["TIME_LIMIT"]))
            # Now get all the results
            for i, job in enumerate(job_list):
                try:
                    job_results.append(job.result())
                except Exception as error:
                    print(error.args, i)
            
            prefix = galprime.gen_out_prefix(current_bin.bin_params)
            index_prefix = galprime.gen_index_prefix(current_bin.bin_params, self.config)

            bare_hdulist = fits.HDUList()
            bgadded_hdulist = fits.HDUList()
            bgsub_hdulist = fits.HDUList()

            # Clean the jobs of any that just so happened to fail and return junk profiles
            cleaned_containers = []
            for container in job_results:
                good = True
                for prof in (container.model_profile, container.bgadded_profile, container.bgsub_profile):
                    if len(prof) == 0:
                        good = False
                        break
                if good:
                    cleaned_containers.append(container)
            print("Profiles cleaned", len(job_results), len(cleaned_containers))
            job_results = cleaned_containers
                    
                        
            for i, container in enumerate(job_results):
                bare_hdulist.append(fits.BinTableHDU(data=galprime.table_from_isolist(container.model_profile)))
                bgadded_hdulist.append(fits.BinTableHDU(data=galprime.table_from_isolist(container.bgadded_profile)))
                bgsub_hdulist.append(fits.BinTableHDU(data=galprime.table_from_isolist(container.bgsub_profile)))

            try:
                bare_hdulist.writeto(bare_dir + index_prefix + "bare.fits", overwrite=True)
                bgadded_hdulist.writeto(bgadded_dir + index_prefix + "bgadded.fits", overwrite=True)
                bgsub_hdulist.writeto(bgsub_dir + index_prefix + "bgsub.fits", overwrite=True)
            except Exception as error:
                print("Error Generating Data Products (Not good!):", error.args)

            try:
                bare_profiles = [galprime.table_from_isolist(container.model_profile) for container in job_results]
                bare_median_table = galprime.boostrap_median(bare_profiles)
                bare_median_table.write(bare_median_dir + index_prefix + "medians.fits", overwrite=True)
            except Exception as error:
                print("Error Bootstrapping Bare:", error.args)
            try:    
                bgadded_profiles = [galprime.table_from_isolist(container.bgadded_profile) for container in job_results]
                bgadded_median_table = galprime.boostrap_median(bgadded_profiles)
                bgadded_median_table.write(bgadded_median_dir + index_prefix + "medians.fits", overwrite=True)
            except Exception as error:
                print("Error Bootstrapping Bgadded:", error.args)
            try:
                bgsub_profiles = [galprime.table_from_isolist(container.bgsub_profile) for container in job_results]
                bgsub_median_table = galprime.boostrap_median(bgsub_profiles)
                bgsub_median_table.write(bgsub_median_dir + index_prefix + "medians.fits", overwrite=True)
            except Exception as error:
                print("Error Bootstrapping Bgsub:", error.args)


def gprime_single(gprime_obj, current_bin, kde):
    metadata = {}
    metadata["NAMES"] = current_bin.object_column_names
    metadata["BIN_PARAMS"] = current_bin.bin_params

    # Select a background and add the background info to the metadata
    bg_index = randint(len(gprime_obj.backgrounds.cutouts))
    metadata["BG_INDEX"] = bg_index
    metadata.update(gprime_obj.backgrounds.cutout_data[bg_index])

    # Get the best PSF is the user has supplied PSFs to the pipeline
    # Otherwise, just use a simple Gaussian2D kernel from astropy
    if gprime_obj.psfs is None:
        psf = Gaussian2DKernel(x_stddev=1.5, y_stddev=1.5)
    else:
        psf = galprime.get_closest_psf(gprime_obj.psfs, metadata["RA"], metadata["DEC"])

    # Generate our container
    container = galprime.GalPrimeContainer(config=gprime_obj.config, kde=kde, psf=psf,
                                    background_cutout=gprime_obj.backgrounds.cutouts[bg_index], 
                                    metadata=metadata)
    container.process_object(plot=True)
    return container
