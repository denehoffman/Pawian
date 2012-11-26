
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <cassert>



#include "TROOT.h"
// #include <TSystem.h>
//PiPiSWaveTMatrix class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

class TFile;
class TH1F;
class TH2F;


class PiPiSWaveTMatrix {

public:

  // create/copy/destroy:

  ///Constructor 
  PiPiSWaveTMatrix();


  /** Destructor */
  virtual ~PiPiSWaveTMatrix();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _invPiPiMassH1;
  TH1F* _invPiPiMassRelH1;
  TH1F* _absT00RelH1;
  TH1F* _absS00RelH1;
  TH1F* _sqrT00RelSigmaPoleH1;
  TH2F* _pipiPhaseSpaceFactorH2;
  TH2F* _pipipipiPhaseSpaceFactorH2;
  TH2F*  _argandH2;
  TH2F*  _argandRelH2;
  TH2F* _phaseShiftH2;
  TH2F* _phaseShiftRelH2;
};


