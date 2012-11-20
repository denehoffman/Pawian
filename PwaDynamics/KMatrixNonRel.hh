// KMatrixNonRel class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixNonRel.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

class KPole;
class AbsPhaseSpace;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixNonRel : public KMatrixBase {

public:

  /// Constructor 
  KMatrixNonRel(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs); 

  /// Destructor
  virtual ~KMatrixNonRel();

  virtual void evalMatrix(const double mass);

protected:
};
//_____________________________________________________________________________


