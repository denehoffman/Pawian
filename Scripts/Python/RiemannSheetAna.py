#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from plotstyle import SetupPlotStyle
SetupPlotStyle()
from scipy.interpolate import griddata

import RiemannSheetAna_py as rs
theRiemannSheetAna = rs.RiemannSheetAna_py()

# result = theRiemannSheetAna.calcTMat(1.3, -0.2)

# theRiemannSheetAna.SetParamValue("a21700Mass", 1.4)

colors=['#840000', '#ffab0f', '#b6c406', '#89a203', '#01889f', '#014182', '#9e0168', 'red']
colorm = LinearSegmentedColormap.from_list('my_map', colors, N=200)

x = []
y = []
z = []

re_Min = 0.8
re_Max = 2.0
re_binning = 0.005
re_nBins = int(round((re_Max-re_Min)/re_binning))

im_Min = -0.01
im_Max = -0.3
im_binning = 0.005
im_nBins = int(round(abs(im_Max-im_Min)/im_binning))

for eReal in np.arange(re_Min, re_Max, re_binning) :
    for eImag in np.arange(im_Min,im_Max, -im_binning) :
        result = theRiemannSheetAna.calcTMat(eReal, eImag)
        x.append(eReal)
        y.append(eImag)
        z.append(result)

xi = np.linspace(re_Min, re_Max, re_nBins)
yi = np.linspace(im_Max, im_Min, im_nBins)
zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')

# fig, ax = plt.subplots(figsize=(10,3))
fig, (ax, ax1) = plt.subplots(2,1)

#ax.set_title("Sheet ----")
ax.set_xlim(re_Min, re_Max)
ax.set_ylim(im_Max, im_Min)

ax.set_xlabel("Re(s)")
ax.set_ylabel("Im(s)")


Norm = cm.colors.Normalize(vmax=max(z), vmin=0.)
surf_plotf = ax.contourf(xi, yi, zi, 10, norm = Norm, cmap=colorm, levels = np.arange(0., max(z), 0.02))
surf_plot = ax.contour(xi, yi, zi, 10, cmap=colorm, linewidths=(0.5), levels = np.arange(0., max(z), 0.02), alpha=0.5)

# cbar = fig.colorbar(surf_plotf, fraction=0.046, pad=0.04)
# cbar.ax.set_ylabel("$|T|^2$")

theRiemannSheetAna.SetParamValue("a21700Mass", 1.4)

x1 = []
y1 = []
z1 = []

for eReal in np.arange(re_Min, re_Max, re_binning) :
    for eImag in np.arange(im_Min,im_Max, -im_binning) :
        result = theRiemannSheetAna.calcTMat(eReal, eImag)
        x1.append(eReal)
        y1.append(eImag)
        z1.append(result)

x1i = np.linspace(re_Min, re_Max, re_nBins)
y1i = np.linspace(im_Max, im_Min, im_nBins)
z1i = griddata((x1, y1), z1, (x1i[None, :], y1i[:, None]), method='cubic')        

ax1.set_xlim(re_Min, re_Max)
ax1.set_ylim(im_Max, im_Min)

ax1.set_xlabel("Re(s)")
ax1.set_ylabel("Im(s)")

Norm1 = cm.colors.Normalize(vmax=max(z1), vmin=0.)
surf_plotf1 = ax1.contourf(x1i, y1i, z1i, 10, norm = Norm1, cmap=colorm, levels = np.arange(0., max(z), 0.02))
surf_plot1 = ax1.contour(x1i, y1i, z1i, 10, cmap=colorm, linewidths=(0.5), levels = np.arange(0., max(z), 0.02), alpha=0.5)

# cbar1 = fig.colorbar(surf_plotf1, fraction=0.046, pad=0.04)
# cbar1.ax1.set_ylabel("$|T|^2$")
plt.show()


# plt.clf()

# theRiemannSheetAna.calcResidue()

# theRiemannSheetAna.calcFVecResidue()
