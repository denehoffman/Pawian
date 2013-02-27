// Flatte class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file Flatte.h
//_____________________________________________________________________________

#include <iostream>
#include <vector>
#include <complex>
#include <utility>
#include <boost/shared_ptr.hpp>

class Particle;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class Flatte {

public:

  /// Constructor 
  Flatte(std::pair<Particle*, Particle*>& decPair1, std::pair<Particle*, Particle*>& decPair2); //decPair1: decay particles
                                                                                                //decPair2: second decay channel    
  Flatte(std::pair <const double, const double>& massPair1, std::pair <const double, const double>& massPair2);
 
  /// Destructor
  ~Flatte();

  complex<double> calcFirstChannel(double currentMass, double mass0, double g1, double g2);
  complex<double> calcSecondChannel(double currentMass, double mass0, double g1, double g2);

protected:
  double _mass11;
  double _mass12;
  double _mass21;
  double _mass22; 
};



