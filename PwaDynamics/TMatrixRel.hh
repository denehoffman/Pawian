// TMatrixRel class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file TMatrixRel.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

class KMatrixBase;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class TMatrixRel : public TMatrixBase {

public:

  /// Constructor 
  TMatrixRel(boost::shared_ptr<KMatrixBase> Kmatrix); 

  /// Destructor
  virtual ~TMatrixRel();

  virtual void evalMatrix(const double mass);

protected:

};
//_____________________________________________________________________________


