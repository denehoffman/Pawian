// AbsPhaseSpace class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file KMatrixBase.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include <complex>
using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class AbsPhaseSpace {

public:

  /// Constructor 
  AbsPhaseSpace(){;} 

  /// Destructor
  virtual ~AbsPhaseSpace(){;}

  // operators:


  // functions:

  virtual complex<double> factor(const double mass)=0;
  virtual complex<double> breakUpMom(const double mass)=0;


protected:

private:

};
//_____________________________________________________________________________


