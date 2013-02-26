// FlatteShape class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file FlatteShape.h
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

class FlatteShape {

public:

  /// Constructor 
  FlatteShape(std::pair<Particle*, Particle*>& decPair1, std::pair<Particle*, Particle*>& decPair2); //decPair1: decay particles
                                                                                                //decPair2: second decay channel    

  /// Destructor
  ~FlatteShape();

  complex<double> calculate(double currentMass, double mass0, double g1, double g2);

protected:
  double _mass11;
  double _mass12;
  double _mass21;
  double _mass22; 
};



