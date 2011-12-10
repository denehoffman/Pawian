#ifndef _JpsiGamKsKlKKHist_H
#define _JpsiGamKsKlKKHist_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"

using JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData;

class JpsiGamKsKlKKEventList;
class AbsJpsiGamKsKlKKLh;
class TFile;
class TH2F;
class TH1F;
class TNtuple;

class JpsiGamKsKlKKHist {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamKsKlKKHist(boost::shared_ptr<const JpsiGamKsKlKKEventList>);
  JpsiGamKsKlKKHist(boost::shared_ptr<AbsJpsiGamKsKlKKLh>, paramGamKsKlKK&);

  /** Destructor */
  virtual ~JpsiGamKsKlKKHist();

  // Getters:
 
protected:


private:
  TFile* _theTFile;
  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;
  
  TH1F* _PhiPhiMassDataHist;
  TH1F* _PhiPhiMassMcHist;
  TH1F* _PhiPhiMassFittedHist;
  
  TH1F* _KpKmMassDataHist;
  TH1F* _KpKmMassMcHist;
  TH1F* _KpKmMassFittedHist;
  TH1F* _KsKlMassDataHist;
  TH1F* _KsKlMassMcHist;
  TH1F* _KsKlMassFittedHist;

  TH1F* _costKs_KsKlHeliDataHist;
  TH1F* _costKs_KsKlHeliMcHist;
  TH1F* _costKs_KsKlHeliFittedHist;
  TH1F* _phiKs_KsKlHeliDataHist;
  TH1F* _phiKs_KsKlHeliMcHist;
  TH1F* _phiKs_KsKlHeliFittedHist;
 
  TH1F*  _costKp_KpKmHeliDataHist; 
  TH1F*  _costKp_KpKmHeliMcHist; 
  TH1F*  _costKp_KpKmHeliFittedHist;
  TH1F*  _phiKp_KpKmHeliDataHist; 
  TH1F*  _phiKp_KpKmHeliMcHist; 
  TH1F*  _phiKp_KpKmHeliFittedHist;
 
  TNtuple* _dataTuple;
  TNtuple* _mcTuple;



  void initRootStuff();
  void plotDalitz(TH2F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight);
  void plotPhiPhi(TH1F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight);
  void plotKsKl(TH1F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight);
  void plotKpKm(TH1F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight);
  void plotCostPhiKs(TH1F* theCostHisto, TH1F* thePhiHisto, const JpsiGamKsKlKKEvtData* theData, double weight);
  void plotCostPhiKp(TH1F* theCostHisto,  TH1F* thePhiHisto, const JpsiGamKsKlKKEvtData* theData, double weight);

  void fillTuple( TNtuple* theTuple, const JpsiGamKsKlKKEvtData* theData, double weight);

};

#endif
