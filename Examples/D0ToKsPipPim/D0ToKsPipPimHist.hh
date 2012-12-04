#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>


#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "qft++/topincludes/tensor.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"

//#include  "Examples/JpsiGamKsKlKK/FitParamErrorMatrix.hh"

class TFile;
class TH2F;
class TH1F;
class TNtuple;
class EvtDataBaseListNew;
//class FitParamErrorMatrix;

class D0ToKsPipPimHist {

public:

  // create/copy/destroy:

  ///Constructor 
  D0ToKsPipPimHist(boost::shared_ptr<const EvtDataBaseListNew>);
  D0ToKsPipPimHist(boost::shared_ptr<AbsLhNew>, fitParamsNew&);
  
  
  /** Destructor */
  virtual ~D0ToKsPipPimHist();

  // Getters:
 
protected:


private:
 
  TFile* _theTFile;
  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;
  
  TH1F* _pippimMassDataHist;
  TH1F* _pippimMassMcHist;
  TH1F* _pippimMassFittedHist;

  TH1F* _KspipMassDataHist;
  TH1F* _KspipMassMcHist;
  TH1F* _KspipMassFittedHist;

  TH1F* _KspimMassDataHist;
  TH1F* _KspimMassMcHist;
  TH1F* _KspimMassFittedHist;

  void initRootStuff();
  void plotDalitz(TH2F* theHisto, EvtDataNew* theData, double weight);
  void plotPiPi(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotKsPip(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotKsPim(TH1F* theHisto, EvtDataNew* theData, double weight);

};

