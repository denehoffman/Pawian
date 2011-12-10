#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKHist.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "Examples/JpsiGamKsKlKK/AbsJpsiGamKsKlKKLh.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKHist::JpsiGamKsKlKKHist(boost::shared_ptr<const JpsiGamKsKlKKEventList> theEvtList) :
  _theTFile(0),
  _PhiPhiMassDataHist(0),
  _PhiPhiMassMcHist(0),
  _PhiPhiMassFittedHist(0),
  _KpKmMassDataHist(0),
  _KpKmMassMcHist(0),
  _KpKmMassFittedHist(0),
  _KsKlMassDataHist(0),
  _KsKlMassMcHist(0),
  _KsKlMassFittedHist(0),
  
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==theEvtList){
    Alert <<"JpsiGamKsKlKKEventList* theEvtList is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  initRootStuff();

  const std::vector<JpsiGamKsKlKKEvtData*> dataList=theEvtList->getDataVecs();

  std::vector<JpsiGamKsKlKKEvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotPhiPhi(_PhiPhiMassDataHist, (*it), 1.  );
      plotKsKl( _KsKlMassDataHist, (*it), 1. );
      plotKpKm( _KpKmMassDataHist, (*it), 1. );
      plotCostPhiKs( _costKs_KsKlHeliDataHist, _phiKs_KsKlHeliDataHist,(*it), 1. );
      plotCostPhiKs( _costKp_KpKmHeliDataHist, _phiKp_KpKmHeliDataHist,(*it), 1. );
      plotCostGam( _costGamCmDataHist,(*it), 1. );      
      fillTuple(_dataTuple, (*it), 1.);

      ++it;
    }

  const std::vector<JpsiGamKsKlKKEvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotPhiPhi(_PhiPhiMassMcHist, (*it), 1.  );
      plotKsKl( _KsKlMassMcHist, (*it), 1. );
      plotKpKm( _KpKmMassMcHist, (*it), 1. );
      plotCostPhiKs( _costKs_KsKlHeliMcHist, _phiKs_KsKlHeliMcHist,(*it), 1. );
      plotCostPhiKp( _costKp_KpKmHeliMcHist, _phiKp_KpKmHeliMcHist, (*it), 1. );
      plotCostGam( _costGamCmMcHist, (*it), 1. );

     fillTuple(_mcTuple, (*it), 1.);

      ++it;
    }
}

JpsiGamKsKlKKHist::JpsiGamKsKlKKHist(boost::shared_ptr<AbsJpsiGamKsKlKKLh> theJpsiGamKsKlKKLh, paramGamKsKlKK& fitParam) :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _PhiPhiMassDataHist(0),
  _PhiPhiMassMcHist(0),
  _PhiPhiMassFittedHist(0),
  _KpKmMassDataHist(0),
  _KpKmMassMcHist(0),
  _KpKmMassFittedHist(0),
  _KsKlMassDataHist(0),
  _KsKlMassMcHist(0),
  _KsKlMassFittedHist(0),
  
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==theJpsiGamKsKlKKLh){
    Alert <<"JpsiGamKsKlKKLh* theJpsiGamKsKlKKLh is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  initRootStuff();

  boost::shared_ptr<const JpsiGamKsKlKKEventList> theEvtList=theJpsiGamKsKlKKLh->getEventList();
  const std::vector<JpsiGamKsKlKKEvtData*> dataList=theEvtList->getDataVecs();

  std::vector<JpsiGamKsKlKKEvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotPhiPhi(_PhiPhiMassDataHist, (*it), 1.  );
      plotKsKl( _KsKlMassDataHist, (*it), 1. );
      plotKpKm( _KpKmMassDataHist, (*it), 1. );
      plotCostPhiKs( _costKs_KsKlHeliDataHist,  _phiKs_KsKlHeliDataHist, (*it), 1. );
      plotCostPhiKp( _costKp_KpKmHeliDataHist, _phiKp_KpKmHeliDataHist, (*it), 1. ); 
      plotCostGam( _costGamCmDataHist, (*it), 1. ); 
     
      fillTuple(_dataTuple, (*it), 1.);
      
      ++it;
    }
  
  const std::vector<JpsiGamKsKlKKEvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotPhiPhi(_PhiPhiMassMcHist, (*it), 1.  );
      plotKsKl( _KsKlMassMcHist, (*it), 1. );
      plotKpKm( _KpKmMassMcHist, (*it), 1. );
      plotCostPhiKs( _costKs_KsKlHeliMcHist, _phiKs_KsKlHeliMcHist, (*it), 1. );
      plotCostPhiKp( _costKp_KpKmHeliMcHist, _phiKp_KpKmHeliMcHist,(*it), 1. );
      plotCostGam( _costGamCmMcHist,(*it), 1. );
      
      double evtWeight= theJpsiGamKsKlKKLh->calcEvtIntensity((*it), fitParam);
      plotDalitz(_dalitzFittedHist, (*it),evtWeight );
      plotPhiPhi(_PhiPhiMassFittedHist, (*it), evtWeight  );
      plotKsKl( _KsKlMassFittedHist, (*it), evtWeight );
      plotKpKm( _KpKmMassFittedHist, (*it), evtWeight );
      plotCostPhiKs( _costKs_KsKlHeliFittedHist, _phiKs_KsKlHeliFittedHist,(*it), evtWeight );
      plotCostPhiKp( _costKp_KpKmHeliFittedHist, _phiKp_KpKmHeliFittedHist,(*it), evtWeight );
      plotCostGam( _costGamCmFittedHist,(*it), evtWeight );
      
      fillTuple(_mcTuple, (*it), evtWeight);
      ++it;
    }
  
  
  double integralData=(double) theEvtList->getDataVecs().size();
  Info <<"No of fit data events  " << integralData ;  // << endmsg;
  
  double integralFitted=(double) theEvtList->getMcVecs().size();
  Info <<"No of fit events " << integralFitted ;  // << endmsg; 

  Info <<"scaling factor  " << integralData/integralFitted ;  // << endmsg; 

  double scaleFactor = integralData/integralFitted;
  
  _dalitzFittedHist->Scale(scaleFactor);
  _PhiPhiMassFittedHist->Scale(scaleFactor);
  _KsKlMassFittedHist->Scale(scaleFactor);
  _KpKmMassFittedHist->Scale(scaleFactor);
  _costKs_KsKlHeliFittedHist->Scale(scaleFactor);
  _phiKs_KsKlHeliFittedHist->Scale(scaleFactor);
  _costKp_KpKmHeliFittedHist->Scale(scaleFactor);
  _phiKp_KpKmHeliFittedHist->Scale(scaleFactor);
  _costGamCmFittedHist->Scale(scaleFactor);  
}

JpsiGamKsKlKKHist::~JpsiGamKsKlKKHist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void JpsiGamKsKlKKHist::initRootStuff()
{
  std::string rootFileName="./JpsiGamKsKlKK.root";

  _theTFile=new TFile(rootFileName.c_str(),"recreate");
  int xbins = 50;
  double xmin=0.8;
  double xmax =10.;
  int ybins=xbins;
  double ymin=xmin;
  double ymax=xmax;
  
  _dalitzDataHist= new TH2F("_dalitzDataHist","Dpl K+K- K+#pi^{0} data",xbins, xmin, xmax, ybins, ymin, ymax );
  _dalitzMcHist= new TH2F("_dalitzMcHist","Dpl K+K- K+#pi^{0} MC",xbins, xmin, xmax, ybins, ymin, ymax);
  _dalitzFittedHist= new TH2F("_dalitzFittedHist","Dpl K+K- K+#pi^{0} fit",xbins, xmin, xmax, ybins, ymin, ymax );

  
  int nbins=50;
   xmin=2;
   xmax=3.2;
  _PhiPhiMassDataHist = new TH1F("_PhiPhiMassDataHist", "_PhiPhiMassDataHist", nbins, xmin, xmax);
  _PhiPhiMassMcHist = new TH1F("_PhiPhiMassMcHist", "_PhiPhiMassMcHist", nbins, xmin, xmax);
  _PhiPhiMassFittedHist = new TH1F("_PhiPhiMassFittedHist", "_PhiPhiMassFittedHist", nbins, xmin, xmax);
  
  xmin=0.8;
  xmax=1.2;
  _KsKlMassDataHist = new TH1F("_KsKlMassDataHist", "_KsKlMassDataHist", nbins, xmin, xmax);
  _KsKlMassMcHist = new TH1F("_KsKlMassMcHist", "_KsKlMassMcHist", nbins, xmin, xmax);
  _KsKlMassFittedHist = new TH1F("_KsKlMassFittedHist", "_KsKlMassFittedHist", nbins, xmin, xmax);
  
  _KpKmMassDataHist = new TH1F("_KpKmMassDataHist", "_KpKmMassDataHist", nbins, xmin, xmax);
  _KpKmMassMcHist = new TH1F("_KpKmMassMcHist", "_KpKmMassMcHist", nbins, xmin, xmax);
  _KpKmMassFittedHist = new TH1F("_KpKmMassFittedHist", "_KpKmMassFittedHist", nbins, xmin, xmax);
  
  _costKs_KsKlHeliDataHist= new TH1F("_costKs_KsKlHeliDataHist", "cos(#Theta_{Ks}) KsKlHeli data", 100, -1., 1.); 
  _costKs_KsKlHeliMcHist= new TH1F("_costKs_KsKlHeliMcHist", "cos(#Theta_{Ks}) KsKlHeli Mc", 100, -1., 1.); 
  _costKs_KsKlHeliFittedHist= new TH1F("_costKs_KsKlHeliFittedHist", "cos(#Theta_{Ks}) KsKlHeli Fit", 100, -1, 1);
  _phiKs_KsKlHeliDataHist= new TH1F("_phiKs_KsKlHeliDataHist", "#Phi_{Ks} KsKlHeli data", 100, -1., 1.); 
  _phiKs_KsKlHeliMcHist= new TH1F("_phiKs_KsKlHeliMcHist", "#Phi_{Ks} KsKlHeli Mc", 100, -1., 1.); 
  _phiKs_KsKlHeliFittedHist= new TH1F("_phiKs_KsKlHeliFittedHist", "#Phi_{Ks} KsKlHeli Fit", 100, -1, 1);

  _costKp_KpKmHeliDataHist= new TH1F("_costKp_KpKmHeliDataHist", "cos(#Theta_{K+}) KsKlHeli data", 100, -1., 1.); 
  _costKp_KpKmHeliMcHist= new TH1F("_costKp_KpKmHeliMcHist", "cos(#Theta_{K+}) KsKlHeli Mc", 100, -1., 1.); 
  _costKp_KpKmHeliFittedHist= new TH1F("_costKp_KpKmHeliFittedHist", "cos(#Theta_{K+}) KsKlHeli Fit", 100, -1, 1);
  _phiKp_KpKmHeliDataHist= new TH1F("_phiKp_KpKmHeliDataHist", "#Phi_{K+} KsKlHeli data", 100, -1., 1.); 
  _phiKp_KpKmHeliMcHist= new TH1F("_phiKp_KpKmHeliMcHist", "#Phi_{K+} KsKlHeli Mc", 100, -1., 1.); 
  _phiKp_KpKmHeliFittedHist= new TH1F("_phiKp_KpKmHeliFittedHist", "#Phi_{K+} KsKlHeli Fit", 100, -1, 1);


  _costGamCmDataHist= new TH1F("_costGamCmDataHist", "cos(#Theta_{#gamma}) CM data", 100, -1., 1.); 
  _costGamCmMcHist= new TH1F("_costGamCmMcHist", "cos(#Theta_{#gamma}) CM Mc", 100, -1., 1.); 
  _costGamCmFittedHist= new TH1F("_costGamCmFittedHist", "cos(#Theta_{#gamma}) CM Fit", 100, -1, 1);
 
   std::string tupleVariables = "mKsKlKpKm:mKsKl:mKpKm:KsKlKpKmCostTheta:gamCosTheta:KsKlCosTheta:KpKmCosTheta:KsCosTheta:KpCosTheta:weight";
   
  
   _dataTuple=new TNtuple("_dataTuple", "data ntuple", tupleVariables.data());
   _mcTuple=new TNtuple("_mcTuple", "mc ntuple", tupleVariables.data());
}

void JpsiGamKsKlKKHist::plotDalitz(TH2F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight)
{
  Vector4<float> V4_KsKlKpKm_HeliPsi;
  Vector4<float> V4_KsKl_HeliPsi = theData->V4_KsKl_HeliPsi  ;
  Vector4<float> V4_KpKm_HeliPsi = theData->V4_KpKm_HeliPsi ;
  Vector4<float> V4_gamma_HeliPsi = theData->V4_gamma_HeliPsi ;
  
  double gphi1 = (V4_gamma_HeliPsi + V4_KsKl_HeliPsi ).M();
  double gphi2 = (V4_gamma_HeliPsi + V4_KpKm_HeliPsi ).M();
  theHisto->Fill( gphi1*gphi1, gphi2*gphi2   ,weight);  
}

void JpsiGamKsKlKKHist::plotPhiPhi(TH1F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight){
  Vector4<float> v4 = theData->V4_KsKlKpKm_HeliPsi;
  theHisto->Fill( v4.M(), weight );
}

void JpsiGamKsKlKKHist::plotKsKl(TH1F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight){
  Vector4<float> v4 = theData->V4_KsKl_HeliPsi; 
  theHisto->Fill( v4.M(), weight ); 
}
void JpsiGamKsKlKKHist::plotKpKm(TH1F* theHisto, const JpsiGamKsKlKKEvtData* theData, double weight){
  Vector4<float> v4 = theData->V4_Kp_HeliKpKm; 
  theHisto->Fill( v4.M(), weight );   
}

void JpsiGamKsKlKKHist::plotCostPhiKs(TH1F* theCostHisto,  TH1F* thePhiHisto, const JpsiGamKsKlKKEvtData* theData, double weight){
  Vector4<float> v4 = theData->V4_Ks_HeliKsKl; 
  theCostHisto->Fill( v4.CosTheta(), weight );
  thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamKsKlKKHist::plotCostPhiKp(TH1F* theCostHisto,  TH1F* thePhiHisto, const JpsiGamKsKlKKEvtData* theData, double weight){
  Vector4<float> v4 = theData->V4_Kp_HeliKpKm; 
  theCostHisto->Fill( v4.CosTheta(), weight );
  thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamKsKlKKHist::plotCostGam(TH1F* theCostHisto, const JpsiGamKsKlKKEvtData* theData, double weight){
  Vector4<float> v4 = theData->V4_gamma_HeliPsi; 
  theCostHisto->Fill( v4.CosTheta(), weight );
}


void  JpsiGamKsKlKKHist::fillTuple( TNtuple* theTuple, const JpsiGamKsKlKKEvtData* theData, double weight)
{
  
  Vector4<float> V4_KsKlKpKm_HeliPsi = theData->V4_KsKlKpKm_HeliPsi  ;
  Vector4<float> V4_KsKl_HeliPsi= theData->V4_KsKl_HeliPsi ;
  Vector4<float> V4_KpKm_HeliPsi= theData->V4_KpKm_HeliPsi ;
  Vector4<float> V4_gamma_HeliPsi= theData->V4_gamma_HeliPsi ;
  
  Vector4<float> V4_KsKl_HeliKsKlKpKm= theData->V4_KsKl_HeliKsKlKpKm ;
  Vector4<float> V4_KpKm_HeliKsKlKpKm= theData->V4_KpKm_HeliKsKlKpKm ;
  Vector4<float> V4_Ks_HeliKsKl= theData->V4_Ks_HeliKsKl ;
  Vector4<float> V4_Kp_HeliKpKm= theData->V4_Kp_HeliKpKm ;
  
  
  theTuple->Fill( V4_KsKlKpKm_HeliPsi.M(),
		  V4_KsKl_HeliPsi.M(),
		  V4_KpKm_HeliPsi.M(),
		  V4_KsKlKpKm_HeliPsi.CosTheta(),
		  V4_gamma_HeliPsi.CosTheta(),
		  V4_KsKl_HeliKsKlKpKm.CosTheta(),
		  V4_KpKm_HeliKsKlKpKm.CosTheta(),
		  V4_Ks_HeliKsKl.CosTheta(),
		  V4_Kp_HeliKpKm.CosTheta(),
		  weight);
}
