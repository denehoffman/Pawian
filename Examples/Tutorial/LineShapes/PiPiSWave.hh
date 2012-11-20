#ifndef _PiPiSWave_H
#define _PiPiSWave_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class TFile;
class TH1F;
class TH2F;
class AbsPhaseSpace;

class PiPiSWave {

public:

  // create/copy/destroy:

  ///Constructor 
  PiPiSWave(std::map< int, std::vector<double> >& gFactors, std::vector<double >& poleMasses, std::vector<  boost::shared_ptr<AbsPhaseSpace> > phpVecs);


  /** Destructor */
  virtual ~PiPiSWave();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _invPiPiMassH1;
  TH1F* _invPiPiMassRelH1;
  TH2F*  _argandH2;
  TH2F*  _argandRelH2;
  TH2F* _phaseShiftH2;
  TH2F* _phaseShiftRelH2;
};

#endif
