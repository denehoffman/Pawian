// FVector class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file FVector.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/PVectorRel.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class FVector : public Matrix< complex<double> > {

public:

  /// Constructor 
  FVector(boost::shared_ptr<KMatrixBase> Kmatrix, boost::shared_ptr<PVectorRel> Pvector);
  FVector(int numRows); 

  /// Destructor
  virtual ~FVector();

  virtual void evalMatrix(const double mass);
  virtual void updateBeta(int i, complex<double> beta) {_Pvector->updateBeta(i, beta);}
  virtual void updateFprod (int i, complex<double> fProd) {;}
  virtual void updateS0prod (double s0prod) {;}
  virtual boost::shared_ptr<KMatrixBase> kMatrix(){return _Kmatrix;}
protected:
  boost::shared_ptr<KMatrixBase> _Kmatrix; 
  boost::shared_ptr<PVectorRel> _Pvector; 
};
//_____________________________________________________________________________


