#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
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

re_Min = 0.8
re_Max = 2.0
re_binning = 0.005
re_nBins = int(round((re_Max-re_Min)/re_binning))

im_Min = -0.01
im_Max = -0.3
im_binning = 0.005
im_nBins = int(round(abs(im_Max-im_Min)/im_binning))

def fill_arrays():
    X = []
    Y = []
    Z = []

    for eReal in np.arange(re_Min, re_Max, re_binning) :
        for eImag in np.arange(im_Min,im_Max, -im_binning) :
            result = theRiemannSheetAna.calcTMat(eReal, eImag)
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
    #surf_plot = _ax.contour(xi, yi, zi, 10, cmap=colorm, linewidths=(0.5), levels = np.arange(0., max(z), 0.02), alpha=0.5)

    # cbar = fig.colorbar(surf_plotf, fraction=0.046, pad=0.04)
    # cbar.ax.set_ylabel("$|T|^2$")

    return

fig, ax = plt.subplots(figsize=(10,3))

fig.subplots_adjust(bottom=0.25)

init_mass = 1.32

ax_mass = fig.add_axes([0.15, 0.1, 0.7, 0.05])
mass_slider = Slider(
        ax_mass,
        label='$M(a_2(1320))$',
        valmin=1.,
        valmax=2.,
        valinit=init_mass
)

def update(val):
    theRiemannSheetAna.SetParamValue("a21320Mass", val)
    x, y, z = fill_arrays()
    plot(x,y,z, ax)
    fig.canvas.draw_idle()
  

mass_slider.on_changed(update)

resetax = fig.add_axes([0.8, 0.01, 0.1, 0.05])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    mass_slider.reset()

button.on_clicked(reset)

plt.show()


# plt.clf()

# theRiemannSheetAna.calcResidue()

# theRiemannSheetAna.calcFVecResidue()
