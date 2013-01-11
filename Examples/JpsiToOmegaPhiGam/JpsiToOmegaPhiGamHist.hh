#ifndef _JpsiToOmegaPhiGamHist_H
#define _JpsiToOmegaPhiGamHist_H

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

//#include  "Examples/JpsiGamPipPimPi0KK/FitParamErrorMatrix.hh"

class TFile;
class TH2F;
class TH1F;
class TNtuple;
//class JpsiToOmegaPhiGamLh;
class EvtDataBaseListNew;
//class FitParamErrorMatrix;

class JpsiToOmegaPhiGamHist {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiToOmegaPhiGamHist(boost::shared_ptr<const EvtDataBaseListNew>);
  JpsiToOmegaPhiGamHist(boost::shared_ptr<AbsLhNew>, fitParamsNew&);
  // JpsiToOmegaPhiGamHist(JpsiGamPipPimPi0KKProdLh* theJpsiGamPipPimPi0KKLh, fitParams& fitParam,  FitParamErrorMatrix* theErrorMatrix );
  // void fill();
  void setMassRange(std::pair<double, double> theMassRange){ _massRange = theMassRange; }
  
  
  /** Destructor */
  virtual ~JpsiToOmegaPhiGamHist();

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
  TH1F* _PipPimPi0MassDataHist;
  TH1F* _PipPimPi0MassMcHist;
  TH1F* _PipPimPi0MassFittedHist;

  TH1F* _costPip_PipPimPi0HeliDataHist;
  TH1F* _costPip_PipPimPi0HeliMcHist;
  TH1F* _costPip_PipPimPi0HeliFittedHist;
  TH1F* _phiPip_PipPimPi0HeliDataHist;
  TH1F* _phiPip_PipPimPi0HeliMcHist;
  TH1F* _phiPip_PipPimPi0HeliFittedHist;
 
  TH1F*  _costKp_KpKmHeliDataHist; 
  TH1F*  _costKp_KpKmHeliMcHist; 
  TH1F*  _costKp_KpKmHeliFittedHist;
  TH1F*  _phiKp_KpKmHeliDataHist; 
  TH1F*  _phiKp_KpKmHeliMcHist; 
  TH1F*  _phiKp_KpKmHeliFittedHist;

  TH1F*  _costGamCmDataHist; 
  TH1F*  _costGamCmMcHist; 
  TH1F*  _costGamCmFittedHist;

  TH1F*  _costPhi_KpKmDataHist; 
  TH1F*  _costPhi_KpKmMcHist; 
  TH1F*  _costPhi_KpKmFittedHist; 

  TH1F*  _phiPhi_KpKmDataHist; 
  TH1F*  _phiPhi_KpKmMcHist; 
  TH1F*  _phiPhi_KpKmFittedHist;

  TH1F*  _chiDataHist; 
  TH1F*  _chiMcHist; 
  TH1F*  _chiFittedHist;
 
  TNtuple* _dataTuple;
  TNtuple* _mcTuple;

  TNtuple* _massIndepTuple;
  std::pair<double, double> _massRange;

  void initRootStuff();
  void plotDalitz(TH2F* theHisto, EvtDataNew* theData, double weight);
  void plotPhiPhi(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotPipPimPi0(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotKpKm(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotCostPhiPip(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  void plotCostPhiKp(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  void plotCostPhi_PhiPhiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, const Vector4<double>& the4Vec, double weight);
  void plotCostGam(TH1F* theCostHisto, EvtDataNew* theData, double weight);

  void plotChi(TH1F* theChiHisto, EvtDataNew* theData, double weight);

  void fillTuple( TNtuple* theTuple, EvtDataNew* theData, double weight);

  double decayAngleChi(const Vector4<double>& v4_p,const Vector4<double>& v4_d1,
		       const Vector4<double>& v4_d2,const Vector4<double>& v4_h1,
		       const Vector4<double>& v4_h2 ) ;


  // JpsiGamPipPimPi0KKProdLh* _theJpsiGamPipPimPi0KKLh;
  // fitParams _fitParam;
  // FitParamErrorMatrix* _errMatrix;

};

#endif
