import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

def SetupPlotStyle():
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{xspace}'
    mpl.rcParams['text.usetex'] = True
    #mpl.rcParams['text.latex.preamble'].append(r'\usepackage{siunitx}')
    #mpl.rcParams['figure.figsize'] = (10, 6)
    #mpl.rc('font', family = 'serif', serif = 'cmr10') 
    mpl.rc('font', family = 'serif') 
    mpl.rcParams['axes.labelsize'] = 10#32
    mpl.rcParams['axes.linewidth'] = 0.5#2
    mpl.rcParams['axes.xmargin'] = 0#2
    mpl.rcParams['figure.figsize'] = (5,5)#(20, 20)
    mpl.rcParams['figure.dpi'] =     200
    mpl.rcParams['grid.alpha'] = 0.5
    mpl.rcParams['font.size'] = 10
    mpl.rcParams['font.weight'] = 400
    mpl.rcParams['legend.frameon'] = True
    mpl.rcParams['legend.handletextpad'] = 0.3
    mpl.rcParams['legend.numpoints'] = 1
    mpl.rcParams['legend.labelspacing'] = 0.15
    mpl.rcParams['legend.fontsize'] = 5#28
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['lines.markeredgewidth'] = 0
    mpl.rcParams['lines.markersize'] = 4#8
    mpl.rcParams['savefig.bbox'] =  'tight'
    mpl.rcParams['savefig.pad_inches'] = 0.1
    mpl.rcParams['xtick.major.size']=  5#14
    mpl.rcParams['xtick.minor.size']=  3.5#7
    mpl.rcParams['xtick.major.width']= 0.5#1.5
    mpl.rcParams['xtick.minor.width']= 0.5#1.5
    mpl.rcParams['xtick.major.pad']=   5#10
    mpl.rcParams['xtick.minor.pad']=   5#10
    mpl.rcParams['xtick.labelsize']=   8#30
    mpl.rcParams['ytick.major.size']=  5#14
    mpl.rcParams['ytick.minor.size']=  3.5#7
    mpl.rcParams['ytick.major.width']= 0.5#1.5
    mpl.rcParams['ytick.minor.width']= 0.5#1.5
    mpl.rcParams['ytick.major.pad']=   5#10
    mpl.rcParams['ytick.minor.pad']=   5#10
    mpl.rcParams['ytick.labelsize']=   8#30
    mpl.rcParams['axes.titlesize']=   6#30

