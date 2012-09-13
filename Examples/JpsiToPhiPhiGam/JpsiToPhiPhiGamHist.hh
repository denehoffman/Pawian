#ifndef _JpsiToPhiPhiGamHist_H
#define _JpsiToPhiPhiGamHist_H

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
//class JpsiToPhiPhiGamLh;
class EvtDataBaseListNew;
//class FitParamErrorMatrix;

class JpsiToPhiPhiGamHist {

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiToPhiPhiGamHist(boost::shared_ptr<const EvtDataBaseListNew>);
  JpsiToPhiPhiGamHist(boost::shared_ptr<AbsLhNew>, fitParamsNew&);
  // JpsiToPhiPhiGamHist(JpsiGamKsKlKKProdLh* theJpsiGamKsKlKKLh, fitParams& fitParam,  FitParamErrorMatrix* theErrorMatrix );
  // void fill();
  void setMassRange(std::pair<double, double> theMassRange){ _massRange = theMassRange; }
  
  
  /** Destructor */
  virtual ~JpsiToPhiPhiGamHist();

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
  void plotKsKl(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotKpKm(TH1F* theHisto, EvtDataNew* theData, double weight);
  void plotCostPhiKs(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  void plotCostPhiKp(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight);
  void plotCostPhi_PhiPhiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, const Vector4<double>& the4Vec, double weight);
  void plotCostGam(TH1F* theCostHisto, EvtDataNew* theData, double weight);

  void plotChi(TH1F* theChiHisto, EvtDataNew* theData, double weight);

  void fillTuple( TNtuple* theTuple, EvtDataNew* theData, double weight);

  double decayAngleChi(const Vector4<double>& v4_p,const Vector4<double>& v4_d1,
		       const Vector4<double>& v4_d2,const Vector4<double>& v4_h1,
		       const Vector4<double>& v4_h2 ) ;


  // JpsiGamKsKlKKProdLh* _theJpsiGamKsKlKKLh;
  // fitParams _fitParam;
  // FitParamErrorMatrix* _errMatrix;

};

#endif
