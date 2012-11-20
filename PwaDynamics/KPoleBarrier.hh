// KPoleBarrier class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KPoleBarrier.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include <iostream>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
using namespace std;

class AbsPhaseSpace;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KPoleBarrier : public KPole {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
//   KPoleBarrier() : Matrix<double>::Matrix() {}

  /// Constructor 
  KPoleBarrier(vector<double>& g_i, double mass_0, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, int orbMom); 

  /// Copy Constructor
  // KPoleBarrier(const KPoleBarrier &theCopy);

  /// Destructor
  virtual ~KPoleBarrier();

  // operators:


  // functions:

  virtual void evalMatrix(const double mass);


  
protected:
  vector<boost::shared_ptr<AbsPhaseSpace> > _phpVecs;
  int _orbMom;
};
//_____________________________________________________________________________
