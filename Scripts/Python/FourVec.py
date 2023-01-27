#!/usr/bin/env python
#import sys

#sys.path.insert(1, '/nfs/data/duldul/bertram/PawianPyWrapper/Pawian/lib/')
#print sys.path

import FourVec_py as fv

Vec = fv.FourVec_py()
Vec.greet()
Vec.set(1.,1.,1.,3.)
print(Vec.mag())
