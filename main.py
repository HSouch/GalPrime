import galprime
import sys
from astropy.table import Table
from numpy.random import normal

config = galprime.load_config_file('../newconfig.galprime')
mag_table = "../COSMOS_V9_MAGS.fits"
# galprime.print_config(config)
# sys.exit()

# Load in backgrounds and psfs
bgs = galprime.cutouts_from_file("../cutouts.fits")
psfs = galprime.cutouts_from_file('../i_psfs.fits')

# TODO change this to the actual COSMOS mag table once I have it available to me
mags = normal(loc=25, scale=1, size=500)
mag_kde = galprime.gaussian_kde(mags)

gp = galprime.GPrime(config='../newconfig.galprime', backgrounds=bgs, mag_kde=mag_kde, psfs=psfs)

# test_container.save_outputs(filename="GPrime_" + str(test_container.metadata["ID"]) + ".fits")

if __name__ ==  '__main__':
    test_container = gp.pipeline(max_bins=1, debug=True, process_method=galprime.gprime_single)

