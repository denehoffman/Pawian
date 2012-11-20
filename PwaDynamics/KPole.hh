// KPole class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KPole.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include <iostream>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
using namespace std;

class AbsPhaseSpace;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KPole : public Matrix< complex<double> > {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
//   KPole() : Matrix<double>::Matrix() {}

  /// Constructor 
  KPole(vector<double>& g_i, double mass_0); 

  /// Copy Constructor
  // KPole(const KPole &theCopy);

  /// Destructor
  virtual ~KPole();

  // operators:


  // functions:

  virtual void evalMatrix(const double mass);


  virtual double poleMass() {return _poleMass;}
  virtual vector<double> gFactors() {return _g_i;}

  void updatePoleMass (double newPoleMass) {_poleMass=newPoleMass;}
  
protected:
  vector<double> _g_i;
  double _poleMass;
};
//_____________________________________________________________________________
