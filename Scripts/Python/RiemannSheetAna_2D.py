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
    alpha=0.6,
    antialiased=True,
    rstride=1,
    cstride=1
)

plot_style2 = dict(
    linewidth=0,
    alpha=1.,
    antialiased=True,
    rstride=1,
    cstride=1
)

fig, ax = plt.subplots(
    figsize=(3, 3),
    subplot_kw={"projection": "3d"},
    #tight_layout=True,
)


colors=['#840000', '#ffab0f', '#b6c406', '#89a203', '#01889f', '#014182', '#9e0168', 'red']
colorm = LinearSegmentedColormap.from_list('my_map', colors, N=200)

re_Min = 0.4
re_Max = 2.0
re_binning = 0.005
re_nBins = int(round((re_Max-re_Min)/re_binning))

#im_Min = -0.01
#im_Max = -0.3
im_Min = 0.3
im_Max = 0.0
im_binning = 0.005
im_nBins = int(round(abs(im_Max-im_Min)/im_binning))

def fill_arrays(reMin, reMax, imMin, imMax):

    X = np.arange(reMin, reMax, re_binning)
    Y = np.arange(imMin,imMax, -im_binning)
    X, Y = np.meshgrid(X, Y)
    Z = np.empty([len(X), len(X[0])])

    for i in range(0, len(X)):
        dummy = []
        for j in range(0, len(X[0])):
            result = theRiemannSheetAna.calcTMat(X[i][j], Y[i][j])
            Z[i][j] = result
    return X, Y, Z

def fill_arrays_line(reMin, reMax):

    X = np.arange(reMin, reMax, re_binning)
    Y = np.arange(-0.005, 0.01, 0.005)
    X, Y = np.meshgrid(X, Y)
    Z = np.empty([len(X), len(X[0])])

    for i in range(0, len(X)):
        dummy = []
        for j in range(0, len(X[0])):
            result = theRiemannSheetAna.calcTMat(X[i][j], Y[i][j])
            Z[i][j] = result
    return X, Y, Z

def plot2D(x, y, z, x2, y2, z2, _ax, **kwargs):

    #_ax.set_xlim(re_Min, re_Max)
    #_ax.set_ylim(im_Max, im_Min)
    #_ax.cla()
    _ax.set_zlim(0., 1.5*max(z2[0]))

    _ax.set_xlabel("Re(s)")
    _ax.set_ylabel("Im(s)")

    _ax.plot_surface(x, y, z, **plot_style, **kwargs)
    _ax.plot_surface(x2, y2, z2, **plot_style2, color="red")
    _ax.view_init(elev=10., azim=60)

    return


init_val = 1.32
init_par = "a21320Mass"


#ax = fig.add_subplot(111, projection='3d')
theRiemannSheetAna.SetParamValue(init_par, init_val)
theRiemannSheetAna.SetSheet("--")
x, y, z = fill_arrays(re_Min, re_Max, im_Min, im_Max)
x2, y2, z2 = fill_arrays_line(re_Min, re_Max)
plot2D(x, y, z, x2, y2, z2, ax, color="#02ccfe")
im_Min = -0.0
im_Max = -0.3
theRiemannSheetAna.SetSheet("++")
x_2, y_2, z_2 = fill_arrays(re_Min, re_Max, im_Min, im_Max)
x2_2, y2_2, z2_2 = fill_arrays_line(re_Min, re_Max)
plot2D(x_2, y_2, z_2, x2_2, y2_2, z2_2, ax, color="#89a203")

ax_mass = fig.add_axes([0.15, 0.1, 0.7, 0.05])
mass_slider = Slider(
        ax_mass,
        label='$M(a_2(1320))$',
        valmin=1.,
        valmax=2.,
        valinit=init_val
)


fig.subplots_adjust(bottom=0.25)

def update(val):
    #plt.close(fig)
    theRiemannSheetAna.SetParamValue("a21320Mass", val)
    ax.cla()
    theRiemannSheetAna.SetSheet("--")
    im_Min = 0.3
    im_Max = 0.
    x, y, z = fill_arrays(re_Min, re_Max, im_Min, im_Max)
    x2, y2, z2 = fill_arrays_line(re_Min, re_Max)
    plot2D(x, y, z, x2, y2, z2, ax, color="#02ccfe")
    im_Min = -0.0
    im_Max = -0.3
    theRiemannSheetAna.SetSheet("++")
    x_2, y_2, z_2 = fill_arrays(re_Min, re_Max, im_Min, im_Max)
    x2_2, y2_2, z2_2 = fill_arrays_line(re_Min, re_Max)
    plot2D(x_2, y_2, z_2, x2_2, y2_2, z2_2, ax, color="#89a203")

#    fig.canvas.draw_idle()


mass_slider.on_changed(update)

resetax = fig.add_axes([0.8, 0.01, 0.1, 0.05])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    mass_slider.reset()

button.on_clicked(reset)
update(init_val)

plt.show()


