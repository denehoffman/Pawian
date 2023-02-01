#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import awkward as ak

import HeliVec_py as hv

theHelVec = hv.HeliVec_py()

theHelVec.setMother( 0.3, 0.4, 0.5, 4. )
theHelVec.setDaughter( 0.05, 0.2, 0.1, 0.5 )
theHelVec.setRef( 0.5, -0.2, -0.1, 2.8 )
theHelVec.setMotherRef( 0., 0., 1.0, 3.5 ) 

res = theHelVec.calcHeliVec()

mother = theHelVec.getMother()
daughter= theHelVec.getDaughter()
motherRef = theHelVec.getMotherRef()
ref = theHelVec.getRef()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

origin = [0,0,0]
X, Y, Z = zip(origin)

ax.quiver(X, Y, Z, res[0]/res[3], res[1]/res[3], res[2]/res[3], color='red')

ax.quiver(X, Y, Z, mother[0]/mother[3], mother[1]/mother[3], mother[2]/mother[3], color='orange')

ax.quiver(X, Y, Z, motherRef[0]/motherRef[3], motherRef[1]/motherRef[3], motherRef[2]/motherRef[3], color='fuchsia')
ax.quiver(X, Y, Z, daughter[0]/daughter[3], daughter[1]/daughter[3], daughter[2]/daughter[3], color='blue')
ax.quiver(X,Y,Z, ref[0]/ref[3], ref[1]/ref[3], ref[2]/ref[3], color='green')

ax.set_xlim([-0.5,0.5])
ax.set_ylim([-0.5,0.5])
ax.set_zlim([-0.5,0.5])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
