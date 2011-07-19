#ifndef _BwShape_H
#define _BwShape_H

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

class BwShape {

public:

  // create/copy/destroy:

  ///Constructor 
  BwShape(double MassRes, double MassWidth, double MassDec1, double MassDec2, unsigned int Lmax, double deltaMass);


  /** Destructor */
  virtual ~BwShape();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  std::map <unsigned int, TH1F* > _histMap;
  std::map <unsigned int, TH2F* > _argandHistMap;
  std::map <unsigned int, TH2F* > _phaseHistMap;
};

#endif
