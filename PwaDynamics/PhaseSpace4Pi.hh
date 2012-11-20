// PhaseSpace4Pi class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 
//_____________________________________________________________________________
// @file PhaseSpace4Pi.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include "PwaDynamics/AbsPhaseSpace.hh"

using namespace std;

//_____________________________________________________________________________
//_____________________________________________________________________________

class PhaseSpace4Pi: public AbsPhaseSpace {

public:

  /// Constructor 
  PhaseSpace4Pi(); 

  /// Destructor
  virtual ~PhaseSpace4Pi();

  // operators:


  // functions:

  virtual complex<double> factor(const double mass);
  virtual complex<double> breakUpMom(const double mass);


protected:

private:
  const double _piMass;
};
//_____________________________________________________________________________


