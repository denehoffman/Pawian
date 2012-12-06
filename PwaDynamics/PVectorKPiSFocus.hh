// PVectorKPiSFocus class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file PVectorKPiSFocus.h
//_____________________________________________________________________________

#include "PwaDynamics/PVectorRel.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "PwaDynamics/KMatrixKPiSFocus.hh"

//#include <boost/multi_array.hpp>

class PPole;
class AbsPhaseSpace;

//typedef boost::multi_array< complex<double>, 2> array_type_2dc;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class PVectorKPiSFocus : public PVectorRel {

public:

  /// Constructor 
  PVectorKPiSFocus(boost::shared_ptr<KMatrixKPiSFocus> kMatrix);

  /// Destructor
  virtual ~PVectorKPiSFocus();

  virtual void evalMatrix(const double mass);
  virtual void updateAprod (int i, double aProd);
  virtual void updateBprod (int i, double bProd);
  virtual void updateCprod (int i, double cProd);
  virtual void updatePhaseprod (int i, double phaseProd);

protected:
  std::vector< double > _aProdVec;
  std::vector< double > _bProdVec;
  std::vector< double > _cProdVec;
  std::vector< double > _phaseProdVec;
  boost::shared_ptr<KMatrixKPiSFocus> _kMatrix;
  void init2IsoSpin1();
  void init2IsoSpin3();
};
//_____________________________________________________________________________


