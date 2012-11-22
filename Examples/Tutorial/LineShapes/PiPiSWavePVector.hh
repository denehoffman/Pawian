// PiPiSWavePVector class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once 

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

class PiPiSWavePVector {

public:

  // create/copy/destroy:

  ///Constructor 
  PiPiSWavePVector();


  /** Destructor */
  virtual ~PiPiSWavePVector();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _absF0H1;
  TH1F* _sqrF0H1;
  TH2F*  _argandH2;
  TH2F* _phaseShiftH2;
};


