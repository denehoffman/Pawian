#ifndef _OmegaPiHist_H
#define _OmegaPiHist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"

using OmegaPiData::OmPiEvtData;

class OmegaPiEventList;
class OmegaPiLh;
class TFile;
class TH2F;
class TH1F;


class OmegaPiHist {

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiHist(boost::shared_ptr<const OmegaPiEventList>);
  OmegaPiHist(boost::shared_ptr<OmegaPiLh>, OmegaPiData::fitParamVal&);

  /** Destructor */
  virtual ~OmegaPiHist();

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
 
  unsigned _jmax;
  unsigned _pbarmom;

  void initRootStuff();
  void plotCosOmegaHeli(TH1F* theHisto, const OmPiEvtData* theEvtData, double weight);
  void plotCosPi0FromOmegaHeli(TH1F* theHisto, const OmPiEvtData* theEvtData, double weight);
  void plotTreimanYang(TH1F* theHisto, const OmPiEvtData* theEvtData, double weight);
};

#endif
