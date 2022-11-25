import matplotlib as mpl

def load_pyplot_style():
    """ Load the style used for galprime plotting 
    Notably has a thicker border and tick marks 
    """
    mpl.rc('text', usetex=True)

    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['xtick.minor.size'] = 4
    mpl.rcParams['xtick.minor.width'] = 1
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['ytick.major.width'] = 3
    mpl.rcParams['ytick.minor.size'] = 4
    mpl.rcParams['ytick.minor.width'] = 1
    mpl.rcParams['axes.linewidth'] = 1.5

    mpl.rc('xtick', labelsize=18)
    mpl.rc('ytick', labelsize=18)

    mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

    font = {'family' : 'serif',
        'weight': 'bold',
        'size': 20}

    mpl.rc('font', **font)



