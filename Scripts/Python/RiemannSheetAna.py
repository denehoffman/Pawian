#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import awkward as ak

import RiemannSheetAna_py as rs

theRiemannSheetAna = rs.RiemannSheetAna_py()

result = theRiemannSheetAna.calcTMat(1.3, -0.2)

print(result)

for eReal in np.arange(0.8,2.0, 0.01) :
    for eImag in np.arange(-0.01,-0.3, -0.01) :
        result = theRiemannSheetAna.calcTMat(eReal, eImag)
        print(eReal, eImag, result)
