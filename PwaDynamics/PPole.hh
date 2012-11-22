// PPole class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PPole.h
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

class PPole : public KPole {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
//   PPole() : Matrix<double>::Matrix() {}

  /// Constructor 
  PPole(complex<double>& beta, vector<double>& g_i, double mass_0); 

  /// Copy Constructor
  // PPole(const PPole &theCopy);

  /// Destructor
  virtual ~PPole();

  // operators:


  // functions:

  virtual void evalMatrix(const double mass);

  void updateBeta (complex<double> beta) {_beta=beta;}

  
protected:
  complex<double> _beta;
};
//_____________________________________________________________________________
