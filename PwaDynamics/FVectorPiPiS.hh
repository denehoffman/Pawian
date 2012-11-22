// FVectorPiPiS class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file FVectorPiPiS.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PVectorRel.hh"
#include "PwaDynamics/PVectorSlowCorRel.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class FVectorPiPiS : public FVector {

public:

  /// Constructor 
  FVectorPiPiS(); 

  /// Destructor
  virtual ~FVectorPiPiS();
  virtual void updateFprod (int i, complex<double> fProd){_pVectorCor->updateFprod (i, fProd);}
  virtual void updateS0prod (double s0prod) {_pVectorCor->updateS0prod(s0prod);}
protected:
  boost::shared_ptr<PVectorSlowCorRel> _pVectorCor;
};
//_____________________________________________________________________________


