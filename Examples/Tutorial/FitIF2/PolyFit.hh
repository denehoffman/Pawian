#ifndef _PolyFit_H
#define _PolyFit_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Examples/Tutorial/FitIF2/MinimizerInterface/MIGeneric.hh"

class TFile;
class TGraph;
class TCanvas;
class TRandom;

class PolyFit : public MIGeneric {

public:

  // create/copy/destroy:

  ///Constructor 
  PolyFit(double p0, double p1, double p2, double p3, double sigma);

  PolyFit(const PolyFit&);

  virtual void setStartPar(int num, double *par, double* min, double* max, double* err);

  virtual double evaluator(); // Standard PAWIAN code for the evaluation

  virtual boost::shared_ptr<MIGeneric> clone(){ 
    return boost::shared_ptr<MIGeneric>(new PolyFit(*this));
  };

  /** Destructor */
  virtual ~PolyFit();


  double controlParameter(const std::vector<double>& minPar);
  void drawGraph(double a, double b, double c, double d);
  // Getters:
  std::map <unsigned int, TGraph* > getGraph() const { return _myGraph;};

  std::vector< double > getX() const { return _xValue;};
  std::vector< double > getY() const { return _yValue;};

  const double getSig() const { return _sigma;};
 
protected:


private:
  TFile* _theTFile;
  std::map <unsigned int, TGraph* > _myGraph;

  std::vector< double > _xValue;
  std::vector< double > _yValue;

  double _sigma;

};

#endif
