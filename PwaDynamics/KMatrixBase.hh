// KMatrixBase class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixBase.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "PwaDynamics/AbsPhaseSpace.hh"

class KPole;


using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixBase : public Matrix< complex<double> > {

public:

  /// Constructor 
  KMatrixBase(vector<boost::shared_ptr<KPole> > Kpoles, vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs); 
  KMatrixBase(vector<boost::shared_ptr<AbsPhaseSpace> > phpVecs, int numCols, int numRows);
  KMatrixBase(int numCols, int numRows); 
  /// Destructor
  virtual ~KMatrixBase();

  virtual void evalMatrix(const double mass) {return;}
  virtual vector<boost::shared_ptr<AbsPhaseSpace> > phaseSpaceVec() {return _phpVecs;}
  virtual vector<boost::shared_ptr<KPole> > kpoles() {return _KPoles;}

protected:
  vector<boost::shared_ptr<KPole> > _KPoles;
  vector<boost::shared_ptr<AbsPhaseSpace> > _phpVecs;
};
//_____________________________________________________________________________


