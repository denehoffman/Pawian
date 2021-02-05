#pragma once 

#ifndef INPUTOMNESMATRIXPIPIS5CHANNEL_H
#define INPUTOMNESMATRIXPIPIS5CHANNEL_H

#include <stdio.h>
#include <iostream>
#include <complex>
#include <utility>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "qft++Extension/AbsComplexFunction.hh"

using std::complex;

class InputOmnesMatrixPiPiS5Channel : public AbsComplexFunction{

public:

/// Constructor 
InputOmnesMatrixPiPiS5Channel(unsigned int i, unsigned int j){  cacheAmps(i,j); }

/// Destructor
virtual ~InputOmnesMatrixPiPiS5Channel(){;}
virtual complex<double> eval(double current_s) {return interpolate(current_s);}
virtual complex<double> eval(complex<double> current_s) {return interpolate(current_s);}
  
protected:
  void cacheAmps(unsigned int i, unsigned int j);

  std::ifstream inInputOmnesMatrixPiPiS5Channel;
  int nCalled;
};

#endif /* INPUTOMNESMATRIXPIPIS5CHANNEL_H */
