import galprime
import numpy as np
import os, time

from matplotlib import pyplot as plt

from astropy.io import fits
from astropy.convolution import Gaussian2DKernel

from numpy.random import randint
from scipy.signal import convolve2d

import galprime

from pebble import ProcessPool
from concurrent.futures import TimeoutError


class GalPrimeError(Exception):
    pass


class GPrime():
    """ Class to manage GalPrime simulations """
    
    def __init__(self, config=None, psfs=None, backgrounds=None, mag_kde=None, 
                 table_names=["ID","IMG", "BG_INDEX", "X", "Y", "RA", "DEC", "MAGS", "R50S", "NS", "ELLIPS", "MASS_MED",
                              "ZPHOT", "sfProb", "I_R50", "PA", "R50_PIX", "BG_MEAN", "BG_MED", "BG_STD",
                              "NSEG_BGA", "CENT_BGA", "N_MASKED_BGA", "P_MASKED_BGA", "NSEG_BGSUB", "CENT_BGSUB", 
                              "N_MASKED_BGSUB", "P_MASKED_BGSUB", "T_ELPSD"]):

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

        self.table_names = table_names
        

    def load_bins(self):
        self.binned_objects = galprime.bin_catalog(self.config)
    
    
    def pipeline(self, max_bins=None, mag_kde=None, process_method=None,
                progress_bar=False, debug=False,):
        t_init = time.time()
        # Bin out catalog based on the catalog in the configuration file
        binned_objects = galprime.bin_catalog(self.config)
        max_bins = len(binned_objects.bins) if max_bins is None else max_bins
        verbose = self.config["VERBOSE"]

        if verbose:
            print("Running GalPRIME, version", galprime.__version__)
            print("-" + str(len(binned_objects.bins)), "bins in total. Running GalPRIME on", max_bins, "bins total")
            if self.backgrounds is not None:
                print("-Backgrounds were supplied:", len(self.backgrounds.cutouts), "in total.")


        # Generate the full filestructure, including containers for the individual profiles, and the medians
        # The directories object contains all the required directories in a neat dict format
        directories = galprime.gen_filestructure(self.config["OUT_DIR"])

        # Run through bins and process using the method
        for i in range(max_bins):
            current_bin = binned_objects.bins[i]
            columns = current_bin.columns()
            index_prefix = galprime.gen_index_prefix(current_bin.bin_params, self.config)

            
            # Generate the object KDE that will be used by this bin, and save it to a pickle file 
            kde = galprime.object_kde(columns)
            galprime.save_pickle(kde, directories["ADDITIONAL"] + index_prefix + "_kde.dat")


            # Set up the job pool
            job_list, job_results = [], []
            if verbose:
                print("Working through", self.config["N_MODELS"], "models for bin:", current_bin.bin_params, 
                      "with", self.config["CORES"], "cores." )
            
            # with ProcessPool(max_workers=self.config["CORES"]) as pool:
            #     for i in range(self.config["N_MODELS"]):
            #         job_list.append(pool.schedule(process_method,
            #                                     args=(self, current_bin, kde),
            #                                     timeout=self.config["TIME_LIMIT"]))
            # # Now get all the results
            # for i, job in enumerate(job_list):
            #     try:
            #         job_results.append(job.result())
            #     except Exception as error:
            #         print(" ", error.args, i)

            def task_done(future):
                try:
                    result = future.result()  # blocks until results are ready
                    job_results.append(result)
                except TimeoutError as error:
                    print("Function took longer than %d seconds" % error.args[1])
                except Exception as error:
                    print("Function raised %s" % error)
                    if debug:
                        print(error.traceback)  # traceback of the function


            with ProcessPool(max_workers=self.config["CORES"]) as pool:
                for index in range(self.config["N_MODELS"]):
                    future = pool.schedule(process_method, args=(self, current_bin, kde),
                                           timeout=self.config["TIME_LIMIT"])
                    future.add_done_callback(task_done)


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
            print(" . Profiles cleaned", len(job_results), len(cleaned_containers))
            job_results = cleaned_containers
            

            # Now save all of the data into nice outputs
            bare_hdulist, bgadded_hdulist, bgsub_hdulist = fits.HDUList(), fits.HDUList(), fits.HDUList()
     
            for i, container in enumerate(job_results):
                bare_hdulist.append(fits.BinTableHDU(data=galprime.table_from_isolist(container.model_profile)))
                bgadded_hdulist.append(fits.BinTableHDU(data=galprime.table_from_isolist(container.bgadded_profile)))
                bgsub_hdulist.append(fits.BinTableHDU(data=galprime.table_from_isolist(container.bgsub_profile)))

            try:
                bare_hdulist.writeto(directories["BARE_PROFILES"] + index_prefix + "bare.fits", overwrite=True)
                bgadded_hdulist.writeto(directories["BGADDED_PROFILES"] + index_prefix + "bgadded.fits", overwrite=True)
                bgsub_hdulist.writeto(directories["BGSUB_PROFILES"] + index_prefix + "bgsub.fits", overwrite=True)
            except Exception as error:
                print(" . Error Generating Data Products (Not good!):", error.args)

            try:
                bare_profiles = [galprime.table_from_isolist(container.model_profile) for container in job_results]
                bare_median_table = galprime.boostrap_median(bare_profiles)
                bare_median_table.write(directories["BARE_MEDIANS"] + index_prefix + "medians.fits", overwrite=True)
            except Exception as error:
                print(" . Error Bootstrapping Bare:", error.args)
            try:    
                bgadded_profiles = [galprime.table_from_isolist(container.bgadded_profile) for container in job_results]
                bgadded_median_table = galprime.boostrap_median(bgadded_profiles)
                bgadded_median_table.write(directories["BGADDED_MEDIANS"] + index_prefix + "medians.fits", overwrite=True)
            except Exception as error:
                print(" . Error Bootstrapping Bgadded:", error.args)
            try:
                bgsub_profiles = [galprime.table_from_isolist(container.bgsub_profile) for container in job_results]
                bgsub_median_table = galprime.boostrap_median(bgsub_profiles)
                bgsub_median_table.write(directories["BGSUB_MEDIANS"] + index_prefix + "medians.fits", overwrite=True)
            except Exception as error:
                print(" . Error Bootstrapping Bgsub:", error.args)

            try:
                output_datatable = galprime.construct_infotable(job_results, self.table_names)
                output_datatable.write(directories["ADDITIONAL"] + index_prefix + "_data.fits", 
                                    format="fits", overwrite=True)
            except Exception as error:
                print(" . Error with saving the output table (compare metadata with table_names):", error.args)

            if verbose:
                print("Finished. Time elapsed:", np.round((time.time() - t_init) / 60, 3), "minutes.")

            # ONLY in debugging, you can return the GalPrimeContainer list for the first bin.
            if debug:
                return job_list

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
    container = galprime.GalPrimeContainer(config=gprime_obj.config, kde=kde, psf=np.copy(psf),
                                    background_cutout=np.copy(gprime_obj.backgrounds.cutouts[bg_index]), 
                                    metadata=metadata)
    container.process_object(plot=False)
    return container
