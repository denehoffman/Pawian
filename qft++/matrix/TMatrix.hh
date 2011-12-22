// TMatrix class definition file. -*- C++ -*-
// Copyright 2010 Bertram Kopf
 
#ifndef _TMatrix_H
#define _TMatrix_H
//_____________________________________________________________________________
// @file TMatrix.h
//_____________________________________________________________________________

#include "qft++/matrix/Matrix.hh"
#include "qft++/matrix/KpoleMatrix.hh"
#include <iostream>
#include <vector>
using namespace std;
//_____________________________________________________________________________
//_____________________________________________________________________________

class TMatrix : public Matrix< complex<double> > {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
  //   TMatrix() : Matrix<double>::Matrix() {}

  /// Constructor 
  TMatrix(const vector<KpoleMatrix>& theKpoles); 

  /// Copy Constructor
  //   TMatrix(const TMatrix &theCopy);

  /// Destructor
  virtual ~TMatrix();

  // operators:

//   /// Assignment operator
//   TMatrix& operator=(const Matrix<double> &__){
//     this->Tensor<double>::operator=(__tensor);
//     return *this;
//   }

//   /// Assignment operator
//   TMatrix& operator=(double __x){
//     this->Matrix<double>::operator=(__x);
//     return *this;
//   }

  // functions:

  void updateMatrix(const double mass);
  void updateMatrixRel(const double mass);
  vector< complex<double> > currentRhoFactors() {return _currentRhoFactors;}

protected:

  vector<KpoleMatrix> _Kpoles;
  vector< complex<double> > _currentRhoFactors;
};
//_____________________________________________________________________________

#endif 
