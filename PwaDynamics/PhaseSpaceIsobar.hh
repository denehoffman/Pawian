// PhaseSpaceIsobar class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpaceIsobar.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include "PwaDynamics/AbsPhaseSpace.hh"

using namespace std;

//_____________________________________________________________________________
//_____________________________________________________________________________

class PhaseSpaceIsobar: public AbsPhaseSpace {

public:

  /// Constructor 
  PhaseSpaceIsobar(double mass1, double mass2); 

  /// Destructor
  virtual ~PhaseSpaceIsobar();

  // operators:


  // functions:

  virtual complex<double> factor(const double mass);
  virtual complex<double> breakUpMom(const double mass);


protected:

private:
  double _mass1;
  double _mass2;
};
//_____________________________________________________________________________


