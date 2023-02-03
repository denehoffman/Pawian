#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import awkward as ak

import HeliVec_py as hv

theHelVec = hv.HeliVec_py()

theHelVec.setMother( -1.3, -1.4, 0.5, 4. )
theHelVec.setDaughter( 0.05, 0.2, 0.1, 0.5 )
theHelVec.setRef( 0.5, -0.2, -0.1, 2.8 )
theHelVec.setMotherRef( 0.8, -0.7, -1.0, 3.5 ) 

daughterHeli = theHelVec.calcHeliVec()

mother = theHelVec.getMother()
daughter= theHelVec.getDaughter()
motherRef = theHelVec.getMotherRef()
ref = theHelVec.getRef()

theHelVec.setDaughter(motherRef[0], motherRef[1], motherRef[2], motherRef[3])
motherRefHeli = theHelVec.calcHeliVec()

theHelVec.setDaughter(ref[0], ref[1], ref[2], ref[3])
refHeli  = theHelVec.calcHeliVec()

#fig = plt.figure()
fig = plt.figure(figsize=plt.figaspect(0.5))

# set up the axes for the first plot
ax1 = fig.add_subplot(1, 2, 1, projection='3d')

origin = [0,0,0]
X, Y, Z = zip(origin)

ax1.quiver(X, Y, Z, daughter[0]/daughter[3], daughter[1]/daughter[3], daughter[2]/daughter[3], color='blue', label='daughter')
ax1.quiver(X, Y, Z, mother[0]/mother[3], mother[1]/mother[3], mother[2]/mother[3], color='orange', label='mother')
ax1.quiver(X,Y,Z, ref[0]/ref[3], ref[1]/ref[3], ref[2]/ref[3], color='green', label='ref')
ax1.quiver(X, Y, Z, motherRef[0]/motherRef[3], motherRef[1]/motherRef[3], motherRef[2]/motherRef[3], color='fuchsia', label='motherRef')

ax1.set_xlim([-0.5,0.5])
ax1.set_ylim([-0.5,0.5])
ax1.set_zlim([-0.5,0.5])
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
plt.title('laboratory frame')
ax1.legend(loc='upper right')

ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.quiver(X, Y, Z, daughterHeli[0]/daughterHeli[3], daughterHeli[1]/daughterHeli[3], daughterHeli[2]/daughterHeli[3], color='blue', label='daughterHeli')
ax2.quiver(X, Y, Z, refHeli[0]/refHeli[3], refHeli[1]/refHeli[3], refHeli[2]/refHeli[3], color='green', label='refHeli')
ax2.quiver(X, Y, Z, motherRefHeli[0]/motherRefHeli[3], motherRefHeli[1]/motherRefHeli[3], motherRefHeli[2]/motherRefHeli[3], color='fuchsia', label='motherRefHeli')

ax2.set_xlim([-0.5,0.5])
ax2.set_ylim([-0.5,0.5])
ax2.set_zlim([-0.5,0.5])
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')
plt.title('helicity system (mother rest frame)')
ax2.legend(loc='upper right')

plt.show()
