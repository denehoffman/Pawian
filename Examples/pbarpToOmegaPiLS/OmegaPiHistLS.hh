#ifndef _OmegaPiHistLS_H
#define _OmegaPiHistLS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"
// #include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"

using OmegaPiDataLS::OmPiEvtDataLS;

class AbsOmegaPiEventListLS;
class AbsOmegaPiLhLS;
class TFile;
class TH2F;
class TH1F;


class OmegaPiHistLS {

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiHistLS(boost::shared_ptr<const AbsOmegaPiEventListLS>, const std::string &thePathToRootFile);
  OmegaPiHistLS(boost::shared_ptr<AbsOmegaPiLhLS>, OmegaPiDataLS::fitParamVal&, const std::string &thePathToRootFile);

  /** Destructor */
  virtual ~OmegaPiHistLS();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH1F* _cosOmegaHeliDataHist;
  TH1F* _cosOmegaHeliMcHist;
  TH1F* _cosOmegaHeliFittedHist;
  TH1F* _cosOmegaAccCorHist;
  TH1F* _cosPi0FromOmegaDataHeli;
  TH1F* _cosPi0FromOmegaMcHeli;
  TH1F* _cosPi0FromOmegaFittedHeli;
  TH1F* _cosPi0FromOmegaAccCorHeli;
  TH1F* _treimanYangDataHist;
  TH1F* _treimanYangMcHist;
  TH1F* _treimanYangFittedHist;
  TH1F* _cosPi0FromOmegaDataHeli1; 
  unsigned _lmax;
  unsigned _pbarmom;

  void initRootStuff(const std::string &thePathToRootFile);
  void plotCosOmegaHeli(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotCosPi0FromOmegaHeli(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotCosPi0FromOmegaHeli1(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
  void plotTreimanYang(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight);
};

#endif
