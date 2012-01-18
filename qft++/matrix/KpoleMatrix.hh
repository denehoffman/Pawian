// KpoleMatrix class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf
 
#ifndef _KpoleMatrix_H
#define _KpoleMatrix_H
//_____________________________________________________________________________
// @file KpoleMatrix.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include <iostream>
#include <vector>
using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KpoleMatrix : public Matrix< complex<double> > {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
//   KpoleMatrix() : Matrix<double>::Matrix() {}

  /// Constructor 
  KpoleMatrix(const vector<double>& g_i, const double mass_0, const vector< pair<double, double> >& mDecPair, int orbMom=0); 

  /// Copy Constructor
  KpoleMatrix(const KpoleMatrix &theCopy);

  /// Destructor
  virtual ~KpoleMatrix();

  // operators:


  // functions:

  void updateMatrix(const double mass);
  void updateMatrixRel(const double mass);


  double poleMass() {return _poleMass;}
  vector<double> gFactors() {return _g_i;}
  vector< complex<double> > rhoFactors() {return _rhoPoleFactors;}
  pair<double, double> decayPair(int i) {return _mDecPair[i];}

protected:
  int _orbMom;
  vector<double> _g_i;
  double _poleMass;
  vector< pair<double, double> > _mDecPair;
  vector< complex<double> > _breakupFactors;
  vector< complex<double> > _rhoPoleFactors;
};
//_____________________________________________________________________________

#endif 
