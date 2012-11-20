// KMatrixSlowAdlerCor class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixSlowAdlerCorRel.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

class KPole;
class AbsPhaseSpace;

typedef boost::multi_array<double, 2> array_type_2d;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixSlowAdlerCorRel : public KMatrixBase {

public:

  /// Constructor 
  KMatrixSlowAdlerCorRel(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, boost::shared_ptr<array_type_2d> fscatPtr, double s0_scat, double sAdler0=-0.15, double sAdler=1.); 

  /// Destructor
  virtual ~KMatrixSlowAdlerCorRel();

  virtual void evalMatrix(const double mass);

protected:
  boost::shared_ptr<array_type_2d> _fScatPtr;
  double _s0Scat; 
  double _sAdler0;
  double _sAdler;
};
//_____________________________________________________________________________


