#ifndef _TwoPoles_H
#define _TwoPoles_H

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
class TH1F;
class TH2F;

class TwoPoles {

public:

  // create/copy/destroy:

  ///Constructor 
  TwoPoles(double MassRes1, double Width1, double MassRes2, double Width2, double deltaMass);


  /** Destructor */
  virtual ~TwoPoles();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _massShapeHist;
  TH2F* _argandHist;
  TH2F* _phaseHist;
  TH1F* _massShapeKmatrHist;
  TH2F* _argandKmatrHist;
  TH2F* _phaseKmatrHist;
  TH1F* _massShapeKmatrHistRel;
  TH2F* _argandKmatrHistRel;
  TH2F* _phaseKmatrHistRel;
  double _massLow;
  double _widthLow;
  double _massHigh;
  double _widthHigh;
};

#endif
