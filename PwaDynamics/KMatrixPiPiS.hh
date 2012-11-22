//KMatrixPiPiS class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf

#pragma once 

//_____________________________________________________________________________
// @file KMatrixPiPiS.h
//_____________________________________________________________________________

#include "PwaDynamics/KMatrixSlowAdlerCorRel.hh"
#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>


typedef boost::multi_array< double, 2> array_type_2d;

using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class KMatrixPiPiS : public KMatrixSlowAdlerCorRel {

public:

  /// Constructor 
 KMatrixPiPiS(); 

  /// Destructor
  virtual ~KMatrixPiPiS();

protected:
  void initASParam1900();
private:
};
//_____________________________________________________________________________


