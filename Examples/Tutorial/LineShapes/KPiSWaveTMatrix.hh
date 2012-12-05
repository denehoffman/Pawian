
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <cassert>



#include "TROOT.h"
// #include <TSystem.h>
//KPiSWaveTMatrix class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

class TFile;
class TH1F;
class TH2F;


class KPiSWaveTMatrix {

public:

  // create/copy/destroy:

  ///Constructor 
  KPiSWaveTMatrix();


  /** Destructor */
  virtual ~KPiSWaveTMatrix();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _KPiAmpRealH1;
  TH1F* _KPiAmpImagH1;
};


