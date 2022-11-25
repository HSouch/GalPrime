import sys
from galprime.isophote_l import Ellipse, EllipseGeometry
from numpy import max, pi, log, unravel_index, argmax, ceil
from photutils import data_properties
from tqdm import tqdm



def isophote_fitting(data, config=None, centre_method='standard'):
    """ Wrapper for photutils.isophote methods

    Generates a table of results from isophote fitting analysis. This uses photutils Isophote procedure, which is
    effectively IRAF's Ellipse() method.
    Iterates over many possible input ellipses to force a higher success rate.

    Args:
        data: The input cutout
        config: The configuration parameters, which contains details like STEP and LINEAR for profile extraction.
        centre_method: Which method to use to determine where to place the centre of the first ellipse.
            'standard': place at the exact centre of the cutout
            'max': place at the maximum pixel in the cutout

    Returns:
        The table of results, or an empty list if not fitted successfully.

    """
    # Set-up failsafe in case of strange infinte loops in photutils
    # warnings.filterwarnings("error")

    fail_count, max_fails = 0, 1000
    linear = False if config is None else config["LINEAR"]
    step = 1. if config is None else config["LINEAR_STEP"]
    verbose = False if config is None else config["VERBOSE"]

    # Get centre of image and cutout halfwidth
    if centre_method == 'standard':
        centre = (data.shape[0]/2, data.shape[1]/2)
    elif centre_method == 'max':
        centre = unravel_index(argmax(data), data.shape)
    else:
        centre = (data.shape[0] / 2, data.shape[1] / 2)

    cutout_halfwidth = max((ceil(data.shape[0] / 2), ceil(data.shape[1] / 2)))

    fitting_list = []

    # First, try obtaining morphological properties from the data and fit using that starting ellipse
    try:
        morph_cat = data_properties(log(data))
        r = 2.0
        pos = (morph_cat.xcentroid, morph_cat.ycentroid)
        a = morph_cat.semimajor_sigma.value * r
        b = morph_cat.semiminor_sigma.value * r
        theta = morph_cat.orientation.value

        geometry = EllipseGeometry(pos[0], pos[1], sma=a, eps=(1 - (b / a)), pa=theta)
        flux = Ellipse(data, geometry)
        fitting_list = flux.fit_image(maxit=100, maxsma=cutout_halfwidth, step=step, linear=linear,
                                      maxrit=cutout_halfwidth / 3)
        if len(fitting_list) > 0:
            return fitting_list

    except KeyboardInterrupt:
        sys.exit(1)
    except (RuntimeError, ValueError, OverflowError, IndexError):
        fail_count += 1
        if fail_count >= max_fails:
            return []

    # If that fails, test a parameter space of starting ellipses
    try:
        for angle in range(0, 180, 45):
            for sma in range(2, 26, 5):
                for eps in (0.3, 0.5, 0.9):
                    geometry = EllipseGeometry(float(centre[0]), float(centre[1]), eps=eps,
                                               sma=sma, pa=angle * pi / 180.)
                    flux = Ellipse(data, geometry)
                    fitting_list = flux.fit_image(maxsma=cutout_halfwidth, step=step, linear=linear,
                                                  maxrit=cutout_halfwidth / 3)
                    if len(fitting_list) > 0:
                        return fitting_list

    except KeyboardInterrupt:
        sys.exit(1)
    except (RuntimeError, ValueError, OverflowError, IndexError):

        # print("RuntimeError or ValueError")
        fail_count += 1
        if fail_count >= max_fails:
            return []
    except IndexError:
        fail_count += 1
        if fail_count >= max_fails:
            return []

    return fitting_list
