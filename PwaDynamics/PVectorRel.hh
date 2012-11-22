// PVectorRel class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PVectorRel.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

class PPole;
class AbsPhaseSpace;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class PVectorRel : public KMatrixBase {

public:

  /// Constructor 
  PVectorRel(vector<boost::shared_ptr<PPole> > Ppoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs); 

  /// Destructor
  virtual ~PVectorRel();

  virtual void evalMatrix(const double mass);

  virtual void updateBeta(int i, complex<double> beta);

protected:
  vector<boost::shared_ptr<PPole> > _Ppoles;
};
//_____________________________________________________________________________


