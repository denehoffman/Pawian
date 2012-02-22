#ifndef  _FitParamErrorMatrix_H
#define  _FitParamErrorMatrix_H


#include <iostream>
#include <fstream>

#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MnUserParameters.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


typedef boost::numeric::ublas::matrix<double> boost_matrix;


class FitParamErrorMatrix {
  
public:
  
  FitParamErrorMatrix( ROOT::Minuit2::MnUserCovariance covMatrix, ROOT::Minuit2::MnUserParameters usrParams  );
  FitParamErrorMatrix( int size){_errMatrix = boost_matrix(size,size);}
  FitParamErrorMatrix(std::vector<double> &theData, int size);
  virtual ~FitParamErrorMatrix();
  void Print(std::ostream& os);
  void Write(std::ostream& os);
  std::vector<double> Data();
  double operator() (int row, int col) {return _errMatrix(row,col);}
  

  
private:
  boost_matrix _errMatrix;
  unsigned int _size;

};



#endif
