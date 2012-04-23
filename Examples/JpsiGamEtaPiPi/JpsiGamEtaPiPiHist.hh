#ifndef _JpsiGamEtaPiPiHist_H
#define _JpsiGamEtaPiPiHist_H

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
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsLh.hh"

#include  "Examples/JpsiGamEtaPiPi/FitParamErrorMatrix.hh"

class TFile;
class TH2F;
class TH1F;
class TNtuple;
class JpsiGamEtaPiPiProdLh;
class EvtDataBaseList;
class FitParamErrorMatrix;

class JpsiGamEtaPiPiHist {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiHist(boost::shared_ptr<const EvtDataBaseList>);
  //JpsiGamEtaPiPiHist(boost::shared_ptr<AbsLh>, fitParams&);
  JpsiGamEtaPiPiHist(JpsiGamEtaPiPiProdLh* theJpsiGamEtaPiPiLh, fitParams& fitParam,  FitParamErrorMatrix* theErrorMatrix );
  void fill();
  void setMassRange(std::pair<double, double> theMassRange){ _massRange = theMassRange; }
  
  
  /** Destructor */
  virtual ~JpsiGamEtaPiPiHist();

  // Getters:
 
protected:


private:
  
  TFile* _theTFile;
/*
  TH2F* _dalitzDataHist;
  TH2F* _dalitzMcHist;
  TH2F* _dalitzFittedHist;
*/
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

  TH1F*  _costGamCmDataHist; 
  TH1F*  _costGamCmMcHist; 
  TH1F*  _costGamCmFittedHist;
/*
  TH1F*  _costPhi_KpKmDataHist; 
  TH1F*  _costPhi_KpKmMcHist; 
  TH1F*  _costPhi_KpKmFittedHist; 

  TH1F*  _phiPhi_KpKmDataHist; 
  TH1F*  _phiPhi_KpKmMcHist; 
  TH1F*  _phiPhi_KpKmFittedHist;
*/
  TH1F*  _chiDataHist; 
  TH1F*  _chiMcHist; 
  TH1F*  _chiFittedHist;
 
  TNtuple* _dataTuple;
  TNtuple* _mcTuple;

  TNtuple* _massIndepTuple;
  std::pair<double, double> _massRange;

  void initRootStuff();
  void plotDalitz(TH2F* theHisto, EvtData* theData, double weight);
  void plotEtaPipPim(TH1F* theHisto, EvtData* theData, double weight);
  void plotEtaPi(TH1F* theHisto, EvtData* theData, double weight);
  void plotPipPim(TH1F* theHisto, EvtData* theData, double weight);
  void plotCostPhiEta(TH1F* theCostHisto, TH1F* thePhiHisto, EvtData* theData, double weight);
  void plotCostPhiPip(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtData* theData, double weight);
  //  void plotCostPhi_PhiPhiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, const Vector4<double>& the4Vec, double weight);
  void plotCostGam(TH1F* theCostHisto, EvtData* theData, double weight);

  void plotChi(TH1F* theChiHisto, EvtData* theData, double weight);

  void fillTuple( TNtuple* theTuple, EvtData* theData, double weight);

  double decayAngleChi(const Vector4<double>& v4_p,const Vector4<double>& v4_d1,
		       const Vector4<double>& v4_d2,const Vector4<double>& v4_h1,
		       const Vector4<double>& v4_h2 ) ;


  JpsiGamEtaPiPiProdLh* _theJpsiGamEtaPiPiLh;
  fitParams _fitParam;
  FitParamErrorMatrix* _errMatrix;

};

#endif
