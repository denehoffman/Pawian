#ifndef _DfuncHist_H
#define _DfuncHist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TH1F;
class TFile;

class DfuncHist {

public:

  // create/copy/destroy:

  ///Constructor 
  DfuncHist(Spin J, Spin lam1, Spin lam2);


  /** Destructor */
  virtual ~DfuncHist();

  // Getters:
  void print(std::ostream& os) const; 

protected:


private:
  Spin _J;
  Spin _lam1;
  Spin _lam2;
  TFile* _theTFile;
  TH1F* _dcostHist;
};

#endif
