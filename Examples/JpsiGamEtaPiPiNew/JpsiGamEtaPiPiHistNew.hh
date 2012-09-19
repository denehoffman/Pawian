#ifndef _JpsiGamEtaPiPiHistNew_H
#define _JpsiGamEtaPiPiHistNew_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>


#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
#include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "qft++/topincludes/tensor.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/AbsLhNew.hh"


class TFile;
class TH2F;
class TH1F;
class TNtuple;
class JpsiGamEtaPiPiProdLhNew;
class EvtDataBaseListNew;
class FitParamErrorMatrix;

class JpsiGamEtaPiPiHistNew {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiHistNew(boost::shared_ptr<const EvtDataBaseListNew>,std::pair<double, double> theMassRange);
  JpsiGamEtaPiPiHistNew(boost::shared_ptr<AbsLhNew>, fitParamsNew&, std::pair<double, double> theMassRange);

  //  void setMassRange(std::pair<double, double> theMassRange){ _massRange = theMassRange; }
  double getFitEvents(){
    return _integralFitted;
  };

  
  /** Destructor */
  virtual ~JpsiGamEtaPiPiHistNew();

  // Getters:
 
protected:


private:
  
  TFile* _theTFile;

  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;

  TH1F* _EtaPiPiMassDataHist;
  TH1F* _EtaPiPiMassMcHist;
  TH1F* _EtaPiPiMassFittedHist;
  
  TH1F* _PipPimMassDataHist;
  TH1F* _PipPimMassMcHist;
  TH1F* _PipPimMassFittedHist;
  TH1F* _EtaPiMassDataHist;
  TH1F* _EtaPiMassMcHist;
  TH1F* _EtaPiMassFittedHist;

  TH1F* _costEta_EtaPipHeliDataHist;
  TH1F* _costEta_EtaPipHeliMcHist;
  TH1F* _costEta_EtaPipHeliFittedHist;
  TH1F*  _phiEta_EtaPipHeliDataHist;
  TH1F*  _phiEta_EtaPipHeliMcHist;
  TH1F*  _phiEta_EtaPipHeliFittedHist;
 
  TH1F*  _costPip_PipPimHeliDataHist;
  TH1F*  _costPip_PipPimHeliMcHist;
  TH1F*  _costPip_PipPimHeliFittedHist;
  TH1F*   _phiPip_PipPimHeliDataHist;
  TH1F*   _phiPip_PipPimHeliMcHist;
  TH1F*   _phiPip_PipPimHeliFittedHist;

  TH1F*  _costPiPi_EtaPiPiHeliDataHist;
  TH1F*  _costPiPi_EtaPiPiHeliMcHist;
  TH1F*  _costPiPi_EtaPiPiHeliFittedHist;
  TH1F*   _phiPiPi_EtaPiPiHeliDataHist;
  TH1F*   _phiPiPi_EtaPiPiHeliMcHist;
  TH1F*   _phiPiPi_EtaPiPiHeliFittedHist;

  TH1F*  _costEtaPi_EtaPiPiHeliDataHist;
  TH1F*  _costEtaPi_EtaPiPiHeliMcHist;
  TH1F*  _costEtaPi_EtaPiPiHeliFittedHist;
  TH1F*   _phiEtaPi_EtaPiPiHeliDataHist;
  TH1F*   _phiEtaPi_EtaPiPiHeliMcHist;
  TH1F*   _phiEtaPi_EtaPiPiHeliFittedHist;

  TH1F*  _costGamCmDataHist; 
  TH1F*  _costGamCmMcHist; 
  TH1F*  _costGamCmFittedHist;

 
  TNtuple* _dataTuple;
  TNtuple* _mcTuple;

  TNtuple* _massIndepTuple;
  std::pair<double, double> _massRange;

  double _integralFitted;

  void initRootStuff();
  void rootlogon();
  void plotDalitz(TH2F* theHisto, EvtDataNew* theData, double weight);
  void plotEtaPipPim(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotEtaPi(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotPipPim(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotCostPhiEta(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  void plotCostPhiPip(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  //  void plotCostPhi_PhiPhiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, const Vector4<double>& the4Vec, double weight);
  void plotCostGam(TH1F* theCostHisto, EvtDataNew* theData, double weight);
  void plotCostPhi_PiPiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  void plotCostPhi_EtaPiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight);

  void fillTuple( TNtuple* theTuple, EvtDataNew* theData, double weight);

  fitParamsNew _fitParam;

};

#endif
