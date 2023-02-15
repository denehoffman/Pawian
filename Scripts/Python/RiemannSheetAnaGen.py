import sys
from matplotlib.colors import LinearSegmentedColormap
from plotstyle import SetupPlotStyle
SetupPlotStyle()
from scipy.interpolate import griddata
import sys
import RiemannSheetAna_py as rs

print ("Number of arguments:", len(sys.argv), "arguments")
print ("sys.argv[0]:", sys.argv[0])
print ("sys.argv[1]:", sys.argv[1])
print ("Argument List:", str(sys.argv))

theRiemannSheetAna1 = rs.RiemannSheetAna_py(list(sys.argv))
