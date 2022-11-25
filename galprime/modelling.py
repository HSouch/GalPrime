""" Scripts for generating galaxy models """

from scipy.stats import gaussian_kde
from scipy.special import gamma, gammainc, gammaincinv, kn
from numpy import pi, exp, mgrid, array
from numpy.random import uniform

from astropy.modeling.models import Sersic2D

import galprime

def object_kde(columns):
    return gaussian_kde(columns)


def i_at_r50(mag, n=2, r_50=2, m_0=27):
    """ Get the intensity at the half-light radius """
    b_n = b(n)
    l_tot = 10 ** ((mag - m_0) / -2.5) * (b_n ** (2 * n))
    denom = (r_50 ** 2) * 2 * pi * n * exp(b_n) * gamma(2 * n)
    i_e = l_tot / denom

    return i_e

def b(n, estimate=False):
    """ Get the b_n normalization constant for the sersic profile. 
    From Graham and Driver.
    """
    if estimate:
        return 2 * n - (1 / 3) + (4 / (405 * n)) + (46 / (25515 * (n ** 2)))
    else:
        return gammaincinv(2 * n, 0.5)


def model_from_kde(kde, config=None, mag_kde=None, names=None, seed=None):
    
    cutout_size = 101 if config is None else config["SIZE"]
    arc_conv = 1 if config is None else config["ARC_CONV"]
    
    attempt = array(kde.resample(1, seed=seed), dtype=float)
    params = galprime.ordered_dict(names, attempt)
    
    if mag_kde is not None:
        params["MAGS"] = mag_kde.resample(1, seed=seed)[0]
    
    r_50_pix = params["R50S"] / arc_conv
    
    i_R50 = i_at_r50(params["MAGS"], n=params["NS"], r_50=r_50_pix, m_0=config["ZEROPOINT"])
    theta = uniform(0, 2*pi)
    
    params["I_R50"] = i_R50
    params["PA"] = theta
    params["R50_PIX"] = params["R50S"] / arc_conv
    
    x,y = mgrid[:cutout_size, :cutout_size]
    # Generate the Sersic model from the parameters
    sersic_model = gen_sersic_model(i_R50=i_R50, r_eff=params["R50S"], n=params["NS"], ellip=params["ELLIPS"], 
                                    theta=theta, x_0=cutout_size / 2, y_0=cutout_size/2, 
                                    arc_conv=arc_conv, shape=(cutout_size, cutout_size))
    
    for n in params:
        params[n] = float(params[n])


    return sersic_model, params
    

def gen_sersic_model(i_R50=10, r_eff=5, n=1, ellip=0.1, theta=0, x_0=0, y_0=0, arc_conv=1, shape=(51,51)):
    try:
        sersic_model = Sersic2D(amplitude=i_R50, r_eff=r_eff / arc_conv, n=n, ellip=ellip, theta=theta, x_0=x_0, y_0=y_0)
        x,y = mgrid[:shape[0], :shape[1]]
        z = sersic_model(x,y)

    except ValueError:
        attempt = [i_R50, r_eff, n, ellip, theta, x_0, y_0, arc_conv, shape]
        raise galprime.GalPrimeError("Failed to generate sersic model with the following parameters:\n" + str(attempt))
    
    return z


def check_sersic_params(i_r50, i_r50_max=200):
    pass
