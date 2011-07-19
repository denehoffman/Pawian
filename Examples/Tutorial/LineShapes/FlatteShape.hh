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

  std::string _pType;

};

#endif
