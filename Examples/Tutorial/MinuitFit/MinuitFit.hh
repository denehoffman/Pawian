#ifndef _MinuitFit_H
#define _MinuitFit_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TGraph;
class TCanvas;
class TRandom;

class MinuitFit {

public:

  // create/copy/destroy:

  ///Constructor 
  MinuitFit(double p0, double p1, double p2, double p3, double sigma);


  /** Destructor */
  virtual ~MinuitFit();


  double calcChiSqr(const std::vector<double>& minPar);
  void drawGraph(double a, double b, double c, double d);
  // Getters:
 
protected:


private:
  TFile* _theTFile;
  std::map <unsigned int, TGraph* > _myGraph;

  std::vector< double > _xValue;
  std::vector< double > _yValue;

  double _sigma;

};

#endif
