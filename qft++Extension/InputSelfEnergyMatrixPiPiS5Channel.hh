#pragma once 

#ifndef INPUTSELFENERGYMATRIXPIPIS5CHANNEL_H
#define INPUTSELFENERGYMATRIXPIPIS5CHANNEL_H

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
#include "Utils/PawianConstants.hh"

using std::complex;

class InputSelfEnergyMatrixPiPiS5Channel : public AbsComplexFunction{

public:

/// Constructor 
InputSelfEnergyMatrixPiPiS5Channel(unsigned int i, unsigned int j) : AbsComplexFunction() 
    { cacheAmps(i,j);  }

/// Destructor
virtual ~InputSelfEnergyMatrixPiPiS5Channel(){;}
virtual complex<double> eval(double current_s);
//virtual complex<double> eval(complex<double> current_s);
  
protected:
  void cacheAmps(unsigned int i, unsigned int j);

  virtual complex<double> interpolate(double current_s);

  std::ifstream inInputSelfEnergyMatrixPiPiS5Channel;
  int nCalled;
  double m_1;
  double m_2;
  bool useAnalyticForm=false;

  double kaellen(double current_s, double m1 , double m2);
  complex<double> SE4455(double current_s, double m1, double m2);
};

#endif /* INPUTSelfEnergyMATRIXPIPIS5CHANNEL_H */
