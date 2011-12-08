#ifndef _FlatteShape_H
#define _FlatteShape_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TH1F;
class TH2F;

class FlatteShape {

public:

  // create/copy/destroy:

  ///Constructor 
  FlatteShape(std::string ptype, double g1, double g2);


  /** Destructor */
  virtual ~FlatteShape();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _histShapeLow;
  TH1F* _histShapeHigh;
  TH1F* _histShapeLowTest;
  TH1F* _histShapeHighTest;
  TH1F* _histShapeLowKmatr;
  TH1F* _histShapeHighKmatr;
  TH2F* _argandKmatrFlatteLowHist;
  TH2F* _argandKmatrLowHist;
  TH2F* _argandKmatrHighHist;

  TH1F* _histShapeLowKmatrRel;
  TH1F* _histShapeHighKmatrRel;
  TH2F* _argandKmatrLowRelHist;
  TH2F* _argandKmatrHighRelHist;


  std::string _pType;

};

#endif
