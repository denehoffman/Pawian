//KMatrixKPiSFocus class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @fileKMatrixKPiSFocus.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixBase.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

class KPole;
class AbsPhaseSpace;

typedef boost::multi_array< double, 2> array_type_2d;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixKPiSFocus : public KMatrixBase {

public:

  /// Constructor 
  KMatrixKPiSFocus(int isospin2=1); 
  /// Destructor
  virtual ~KMatrixKPiSFocus();

  virtual void evalMatrix(const double mass);
  virtual double sNorm() {return _sNorm;}

protected:
  boost::shared_ptr<array_type_2d> _aScatPtr;
  boost::shared_ptr<array_type_2d> _bScatPtr;
  boost::shared_ptr<array_type_2d> _cScatPtr;
  double _sAdler0;
  const double _sNorm;

  void init2IsoSpin1();
  void init2IsoSpin3();
};
//_____________________________________________________________________________


