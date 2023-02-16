import sys
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from plotstyle import SetupPlotStyle
SetupPlotStyle()
from scipy.interpolate import griddata
import numpy as np
import RiemannSheetAna_py as rs

print ("Number of arguments:", len(sys.argv), "arguments")
print ("sys.argv[0]:", sys.argv[0])
print ("sys.argv[1]:", sys.argv[1])
print ("Argument List:", str(sys.argv))

colors=['#840000', '#ffab0f', '#b6c406', '#89a203', '#01889f', '#014182', '#9e0168', 'red']
colorm = LinearSegmentedColormap.from_list('my_map', colors, N=200)

re_Min = 0.8
re_Max = 2.0
re_binning = 0.005
re_nBins = int(round((re_Max-re_Min)/re_binning))

im_Min = -0.01
im_Max = -0.3
im_binning = 0.005
im_nBins = int(round(abs(im_Max-im_Min)/im_binning))

def fill_arrays(i):
    X = []
    Y = []
    Z = []

    for eReal in np.arange(re_Min, re_Max, re_binning) :
        for eImag in np.arange(im_Min,im_Max, -im_binning) :
            result = theRiemannSheetAna[i].calcTMat(eReal, eImag)
            X.append(eReal)
            Y.append(eImag)
            Z.append(result)
    return X, Y, Z

def plot(x, y, z, _ax):
    xi = np.linspace(re_Min, re_Max, re_nBins)
    yi = np.linspace(im_Max, im_Min, im_nBins)
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')
    
    ax.set_xlim(re_Min, re_Max)
    ax.set_ylim(im_Max, im_Min)

    ax.set_xlabel("Re(s)")
    ax.set_ylabel("Im(s)")

    Norm = cm.colors.Normalize(vmax=max(z), vmin=0.)
    surf_plotf = _ax.contourf(xi, yi, zi, 10, norm = Norm, cmap=colorm, levels = np.arange(0., max(z), 0.02))
    surf_plot = _ax.contour(xi, yi, zi, 10, cmap=colorm, linewidths=(0.5), levels = np.arange(0., max(z), 0.02), alpha=0.5)

    # cbar = fig.colorbar(surf_plotf, fraction=0.046, pad=0.04)
    # cbar.ax.set_ylabel("$|T|^2$")

    return


theRiemannSheetAna = []

theRiemannSheetAna.append(rs.RiemannSheetAna_py(list(sys.argv)))

sheetSize = theRiemannSheetAna[0].sheetSize()

print ("sheet size:", sheetSize)

def binary(d):
    numbers = []
    i = 1
    for x in range(d):
        numbers.append(i)
        i *= 2
        if i > d:
            i /= 2
            i = int(i)
            break
    numbers.reverse()
    theList = []
    theListstr = ""
    i = 0
    while len(theList) < len(numbers):
        if numbers[i] <= d:
            theList.append("+")
            d -= numbers[i]
        else:
            theList.append("-")
        i += 1
    for x in theList:
        theListstr += str(x)
    return theListstr

sheetList = []

for i in range(0, sheetSize*sheetSize):
    currentSheet=binary(i)
    while len(currentSheet) < sheetSize:
        currentSheet = "-"+ currentSheet
    sheetList.append(currentSheet)
    
sheetList.reverse()
for sheetIdx in range(0,sheetSize*sheetSize):
    print("sheetList[",sheetIdx,"]= ", sheetList[sheetIdx])

for i in range(0, sheetSize*sheetSize):
    if i > 0:
        theRiemannSheetAna.append(rs.RiemannSheetAna_py(list(sys.argv)))
    theRiemannSheetAna[i].SetSheet(sheetList[i])

fig, axes = plt.subplots(nrows=sheetSize, ncols=sheetSize, figsize=(12, 8))

idx = 0
for row in axes:
    for ax in row:
        x, y, z = fill_arrays(idx)
        plot(x, y, z, ax)
        ax.set_title(sheetList[sheetIdx])
        print("fill histo ", idx)
        idx=idx+1

plt.show()
