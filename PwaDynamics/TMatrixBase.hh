// TMatrixBase class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file TMatrixBase.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>


using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class TMatrixBase : public Matrix< complex<double> > {

public:

  /// Constructor 
  TMatrixBase(boost::shared_ptr<KMatrixBase> Kmatrix); 

  /// Destructor
  virtual ~TMatrixBase();

  virtual void evalMatrix(const double mass);

protected:
  boost::shared_ptr<KMatrixBase> _Kmatrix; 

};
//_____________________________________________________________________________


