
import argparse
import galprime
from astropy.table import Table

print("Running with GalPrime: version", galprime.__version__)

# Load in arguments
parser = argparse.ArgumentParser(description="Manager for GalPrime Simulations")

parser.add_argument("CONFIG", type=str, 
                    help="Galprime config file (dump with galprime.dump_default_config_file()).")
parser.add_argument("--mag_table", type=str, help="Optional magnitude table.")
parser.add_argument("--mag_key", type=str, help="Column name for magnitudes")
parser.add_argument("--bgs", type=str, help="Optional supplied backgrounds for GalPrime")
parser.add_argument("--psfs", type=str, help="Supplied PSFs for GalPrime")
parser.add_argument("--max_bins", type=int, help="Maximum bins to run (useful in testing).")
parser.add_argument("--outdir", type=str, help="Output directory (will use one supplied in config file otherwise.")
args = parser.parse_args()
args_dict = args.__dict__


# Load in the configuration file
config_fn = args_dict["CONFIG"]
config = galprime.load_config_file(config_fn)

if args_dict["outdir"] is not None:
    config["OUT_DIR"] = args_dict["outdir"] 

print("All outputs will be saved to", config["OUT_DIR"])


# Load in the optional mags table
if args_dict["mag_table"] is not None:
    mag_table = Table.read(args_dict["mag_table"], format="fits")
    mag_key = args_dict["mag_key"] if args_dict["mag_key"] is not None else "i"

    mags = mag_table[mag_key][mag_table[mag_key] > 0]
    mag_kde = galprime.gaussian_kde(mags)
else:
    mag_kde = None

# Set up backgrounds and psfs
if args_dict["bgs"] is not None:
    bgs = galprime.cutouts_from_file(args_dict["bgs"])
    print(len(bgs.cutouts), "backgrounds loaded.")
else:
    bgs = None

# Load in psfs, either optionally provided first, and if not from the config directly
if args_dict["psfs"] is not None:
    psfs = galprime.cutouts_from_file(args_dict["psfs"])
    print(len(psfs.cutouts), "psfs loaded.")
else:
    psfs  = galprime.cutouts_from_file(config["PSF_FILENAME"])


# Set up max bins
if args_dict["max_bins"] is not None:
    max_bins = int(args_dict["max_bins"])
else:
    max_bins = None    

print()

gp = galprime.GPrime(config=config, backgrounds=bgs, mag_kde=mag_kde, psfs=psfs)
gp.pipeline(max_bins=max_bins, debug=False, process_method=galprime.gprime_single)

