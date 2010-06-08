#ifndef _EtacToapi0Hist_H
#define _EtacToapi0Hist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>
#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Data.hh"

class EtacToapi0EventList;
class TFile;
class TH2F;
class TH1F;


class EtacToapi0Hist {

public:

  // create/copy/destroy:

  ///Constructor 
  EtacToapi0Hist(const EtacToapi0EventList*, const fitParamVal& theParamVal);


  /** Destructor */
  virtual ~EtacToapi0Hist();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;
  TH1F* _invpietaDataHist;
  TH1F* _invpietaMcHist;
  TH1F* _invpietaFittedHist;
  TH1F* _invpipiDataHist;
  TH1F* _invpipiMcHist;
  TH1F* _invpipiFittedHist;

  void initRootStuff();
  void plotDalitz (TH2F* theHisto, const evt4Vec& theEvt4Vecs, double weight);
  void plotInvPiEta (TH1F* theHisto, const evt4Vec& theEvt4Vecs, double weight);
  void plotInvPiPi (TH1F* theHisto, const evt4Vec& theEvt4Vecs, double weight);
};

#endif
