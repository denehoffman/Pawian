#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
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

plot_style = dict(
    linewidth=0,
    alpha=0.9,
    antialiased=True,
    rstride=1,
    cstride=1,
)


#fig, ax = plt.subplots(
#    figsize=(10, 6),
#    subplot_kw={"projection": "3d"},
#    tight_layout=True,
#)


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

def fill_arrays():

    X = np.arange(re_Min, re_Max, re_binning)
    Y = np.arange(im_Min,im_Max, -im_binning)
    X, Y = np.meshgrid(X, Y)
    Z = np.empty([len(X), len(X[0])])

    for i in range(0, len(X)):
        dummy = []
        for j in range(0, len(X[0])):
            result = theRiemannSheetAna.calcTMat(X[i][j], Y[i][j])
            Z[i][j] = result
    return X, Y, Z

def fill_arrays_line():

    X = np.arange(re_Min, re_Max, re_binning)
    Y = np.arange(-0.01, 0., 0.005)
    X, Y = np.meshgrid(X, Y)
    Z = np.empty([len(X), len(X[0])])

    for i in range(0, len(X)):
        dummy = []
        for j in range(0, len(X[0])):
            result = theRiemannSheetAna.calcTMat(X[i][j], Y[i][j])
            Z[i][j] = result
    return X, Y, Z

def plot2D(x, y, z, x2, y2, z2, _ax):

    _ax.set_xlim(re_Min, re_Max)
    _ax.set_ylim(im_Max, im_Min)
    _ax.set_zlim(0., 1.5*max(z2[0]))

    _ax.set_xlabel("Re(s)")
    _ax.set_ylabel("Im(s)")

    _ax.plot_surface(x, y, z, **plot_style, color="#02ccfe")
    _ax.plot_surface(x2, y2, z2, **plot_style, color="red")

    return

#fig, ax = plt.subplots(figsize=(10,3))

#fig.subplots_adjust(bottom=0.25)

fig, ax = plt.subplots(
    figsize=(3, 3),
    subplot_kw={"projection": "3d"},
    tight_layout=True,
)

#ax = fig.add_subplot(111, projection='3d')

init_val = 1.32
init_par = "a21320Mass"

theRiemannSheetAna.SetParamValue(init_par, init_val)
x, y, z = fill_arrays()
x2, y2, z2 = fill_arrays_line()

plot2D(x, y, z, x2, y2, z2, ax)
  

plt.show()


