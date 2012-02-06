#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKHist.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "PwaUtils/KinUtils.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKHist::JpsiGamKsKlKKHist(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
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
  _costPhi_KpKmDataHist(0), 
  _costPhi_KpKmMcHist(0), 
  _costPhi_KpKmFittedHist(0),
  _phiPhi_KpKmDataHist(0), 
  _phiPhi_KpKmMcHist(0), 
  _phiPhi_KpKmFittedHist(0), 
  _chiDataHist(0), 
  _chiMcHist(0), 
  _chiFittedHist(0),
   _dataTuple(0),
  _mcTuple(0)
{
  if(0==theEvtList){
    Alert <<"JpsiGamKsKlKKEventList* theEvtList is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  initRootStuff();

  const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotPhiPhi(_PhiPhiMassDataHist, (*it), 1.  );
      plotKsKl( _KsKlMassDataHist, (*it), 1. );
      plotKpKm( _KpKmMassDataHist, (*it), 1. );
      plotCostPhiKs( _costKs_KsKlHeliDataHist, _phiKs_KsKlHeliDataHist,(*it), 1. );
      plotCostPhiKp( _costKp_KpKmHeliDataHist, _phiKp_KpKmHeliDataHist,(*it), 1. );
      plotCostGam( _costGamCmDataHist,(*it), 1. );
      plotCostPhi_PhiPhiHeli(_costPhi_KpKmDataHist, _phiPhi_KpKmDataHist, (*it)->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm], 1.);        
      fillTuple(_dataTuple, (*it), 1.);
      
      ++it;
    }

  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
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
      plotCostPhi_PhiPhiHeli(_costPhi_KpKmMcHist, _phiPhi_KpKmMcHist, (*it)->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm], 1.);  
     fillTuple(_mcTuple, (*it), 1.);

      ++it;
    }
}

JpsiGamKsKlKKHist::JpsiGamKsKlKKHist(boost::shared_ptr<AbsLh> theJpsiGamKsKlKKLh, fitParams& fitParam) :
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
  _costPhi_KpKmDataHist(0), 
  _costPhi_KpKmMcHist(0), 
  _costPhi_KpKmFittedHist(0),
  _phiPhi_KpKmDataHist(0), 
  _phiPhi_KpKmMcHist(0), 
  _phiPhi_KpKmFittedHist(0), 
  _chiDataHist(0), 
  _chiMcHist(0), 
  _chiFittedHist(0),  
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==theJpsiGamKsKlKKLh){
    Alert <<"JpsiGamKsKlKKLh* theJpsiGamKsKlKKLh is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }

  initRootStuff();

  boost::shared_ptr<const EvtDataBaseList> theEvtList=theJpsiGamKsKlKKLh->getEventList();
  const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotPhiPhi(_PhiPhiMassDataHist, (*it), 1.  );
      plotKsKl( _KsKlMassDataHist, (*it), 1. );
      plotKpKm( _KpKmMassDataHist, (*it), 1. );
      plotCostPhiKs( _costKs_KsKlHeliDataHist,  _phiKs_KsKlHeliDataHist, (*it), 1. );
      plotCostPhiKp( _costKp_KpKmHeliDataHist, _phiKp_KpKmHeliDataHist, (*it), 1. );
      plotCostGam( _costGamCmDataHist, (*it), 1. );
      plotCostPhi_PhiPhiHeli(_costPhi_KpKmDataHist, _phiPhi_KpKmDataHist, (*it)->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm], 1.);  
       plotChi(_chiDataHist, (*it), 1. );      
      fillTuple(_dataTuple, (*it), 1.);
      
      ++it;
    }
  
  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
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
      plotCostPhi_PhiPhiHeli(_costPhi_KpKmMcHist, _phiPhi_KpKmMcHist, (*it)->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm], 1.);  
      plotChi(_chiMcHist, (*it), 1. );   
      
      double evtWeight= theJpsiGamKsKlKKLh->calcEvtIntensity((*it), fitParam);
      plotDalitz(_dalitzFittedHist, (*it),evtWeight );
      plotPhiPhi(_PhiPhiMassFittedHist, (*it), evtWeight  );
      plotKsKl( _KsKlMassFittedHist, (*it), evtWeight );
      plotKpKm( _KpKmMassFittedHist, (*it), evtWeight );
      plotCostPhiKs( _costKs_KsKlHeliFittedHist, _phiKs_KsKlHeliFittedHist,(*it), evtWeight );
      plotCostPhiKp( _costKp_KpKmHeliFittedHist, _phiKp_KpKmHeliFittedHist,(*it), evtWeight );
      plotCostGam( _costGamCmFittedHist,(*it), evtWeight );
      plotCostPhi_PhiPhiHeli(_costPhi_KpKmFittedHist, _phiPhi_KpKmFittedHist, (*it)->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm], evtWeight); 
      plotChi(_chiFittedHist, (*it), evtWeight  );
        
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
  _costPhi_KpKmFittedHist->Scale(scaleFactor);
  _phiPhi_KpKmFittedHist->Scale(scaleFactor);
  _chiFittedHist->Scale(scaleFactor);  
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
  _phiKs_KsKlHeliDataHist= new TH1F("_phiKs_KsKlHeliDataHist", "#Phi_{Ks} KsKlHeli data", 100, -TMath::Pi(), TMath::Pi()); 
  _phiKs_KsKlHeliMcHist= new TH1F("_phiKs_KsKlHeliMcHist", "#Phi_{Ks} KsKlHeli Mc", 100, -TMath::Pi(), TMath::Pi()); 
  _phiKs_KsKlHeliFittedHist= new TH1F("_phiKs_KsKlHeliFittedHist", "#Phi_{Ks} KsKlHeli Fit", 100, -TMath::Pi(), TMath::Pi());

  _costKp_KpKmHeliDataHist= new TH1F("_costKp_KpKmHeliDataHist", "cos(#Theta_{K+}) KsKlHeli data", 100, -1., 1.); 
  _costKp_KpKmHeliMcHist= new TH1F("_costKp_KpKmHeliMcHist", "cos(#Theta_{K+}) KsKlHeli Mc", 100, -1., 1.); 
  _costKp_KpKmHeliFittedHist= new TH1F("_costKp_KpKmHeliFittedHist", "cos(#Theta_{K+}) KsKlHeli Fit", 100, -1, 1);
  _phiKp_KpKmHeliDataHist= new TH1F("_phiKp_KpKmHeliDataHist", "#Phi_{K+} KsKlHeli data", 100, -TMath::Pi(), TMath::Pi()); 
  _phiKp_KpKmHeliMcHist= new TH1F("_phiKp_KpKmHeliMcHist", "#Phi_{K+} KsKlHeli Mc", 100, -TMath::Pi(), TMath::Pi()); 
  _phiKp_KpKmHeliFittedHist= new TH1F("_phiKp_KpKmHeliFittedHist", "#Phi_{K+} KsKlHeli Fit", 100, -TMath::Pi(), TMath::Pi());


  _costGamCmDataHist= new TH1F("_costGamCmDataHist", "cos(#Theta_{#gamma}) CM data", 100, -1., 1.); 
  _costGamCmMcHist= new TH1F("_costGamCmMcHist", "cos(#Theta_{#gamma}) CM Mc", 100, -1., 1.); 
  _costGamCmFittedHist= new TH1F("_costGamCmFittedHist", "cos(#Theta_{#gamma}) CM Fit", 100, -1, 1);

  _costPhi_KpKmDataHist= new TH1F("_costPhi_KpKmDataHist", "cos(#Theta_{#phi}) K+ K- data", 100, -1., 1.);  
  _costPhi_KpKmMcHist= new TH1F("_costPhi_KpKmMcHist", "cos(#Theta_{#phi}) K+ K- Mc", 100, -1., 1.);   
  _costPhi_KpKmFittedHist= new TH1F("_costPhi_KpKmFittedHist", "cos(#Theta_{#phi}) K+ K- Fit", 100, -1., 1.);   

  _phiPhi_KpKmDataHist= new TH1F("_phiPhi_KpKmDataHist", "cos(#Phi_{#phi}) K+ K- data", 100, -TMath::Pi(), TMath::Pi());  
  _phiPhi_KpKmMcHist= new TH1F("_phiPhi_KpKmMcHist", "cos(#Phi_{#phi}) K+ K- Mc", 100, -TMath::Pi(), TMath::Pi());  
  _phiPhi_KpKmFittedHist= new TH1F("_phiPhi_KpKmFittedHist", "cos(#Phi_{#phi}) K+ K- Fit", 100, -TMath::Pi(), TMath::Pi());

  _chiDataHist= new TH1F("_chiDataHist", "#chi data", 90, 0., 90.);
  _chiMcHist= new TH1F("_chiMcHist", "#chi Mc", 90, 0., 90.); 
  _chiFittedHist= new TH1F("_chiFittedHist", "#chi Fit", 90, 0., 90.);  

   std::string tupleVariables = "mKsKlKpKm:mKsKl:mKpKm:KsKlKpKmCostTheta:gamCosTheta:KsKlCosTheta:KpKmCosTheta:KsCosTheta:KpCosTheta:decPlaneChi:testHeli:weight";
   
  
   _dataTuple=new TNtuple("_dataTuple", "data ntuple", tupleVariables.data());
   _mcTuple=new TNtuple("_mcTuple", "mc ntuple", tupleVariables.data());
}

void JpsiGamKsKlKKHist::plotDalitz(TH2F* theHisto,  EvtData* theData, double weight)
{
  Vector4<double>& V4_KsKl_HeliPsi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKl_HeliPsi]  ;
  Vector4<double>& V4_KpKm_HeliPsi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliPsi] ;
  Vector4<double>& V4_gamma_HeliPsi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_gamma_HeliPsi] ;
  
  double gphi1 = (V4_gamma_HeliPsi + V4_KsKl_HeliPsi ).M();
  double gphi2 = (V4_gamma_HeliPsi + V4_KpKm_HeliPsi ).M();
  theHisto->Fill( gphi1*gphi1, gphi2*gphi2   ,weight);  
}

void JpsiGamKsKlKKHist::plotPhiPhi(TH1F* theHisto, EvtData* theData, double weight){
  Vector4<double>& v4 = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  theHisto->Fill( v4.M(), weight );
}

void JpsiGamKsKlKKHist::plotKsKl(TH1F* theHisto, EvtData* theData, double weight){
  Vector4<double>& v4 = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKl_HeliPsi]; 
  theHisto->Fill( v4.M(), weight ); 
}
void JpsiGamKsKlKKHist::plotKpKm(TH1F* theHisto, EvtData* theData, double weight){
  Vector4<double>& v4 = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm]; 
  theHisto->Fill( v4.M(), weight );   
}

void JpsiGamKsKlKKHist::plotCostPhiKs(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtData* theData, double weight){
  Vector4<double>& v4 = theData->FourVecs[enumJpsiGamKsKlKKData::V4_Ks_HeliKsKl]; 
  theCostHisto->Fill( v4.CosTheta(), weight );
  thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamKsKlKKHist::plotCostPhiKp(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtData* theData, double weight){
  Vector4<double>& v4 = theData->FourVecs[enumJpsiGamKsKlKKData::V4_Kp_HeliKpKm]; 
  theCostHisto->Fill( v4.CosTheta(), weight );
  thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamKsKlKKHist::plotCostGam(TH1F* theCostHisto, EvtData* theData, double weight){
  Vector4<double>& v4 = theData->FourVecs[enumJpsiGamKsKlKKData::V4_gamma_HeliPsi]; 
  theCostHisto->Fill( v4.CosTheta(), weight );
}

void JpsiGamKsKlKKHist::plotCostPhi_PhiPhiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, const Vector4<double>& the4Vec, double weight){
  theCostHisto->Fill( the4Vec.CosTheta(), weight);
  thePhiHisto->Fill( the4Vec.Phi(), weight);
}

void JpsiGamKsKlKKHist::plotChi(TH1F* theChiHisto, EvtData* theData, double weight){
//   Vector4<double>& V4_KsKlKpKm_HeliPsi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi]  ;
//   Vector4<double>& V4_Ks_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Ks_HeliPsi] ;
//   Vector4<double>& V4_Kl_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Kl_HeliPsi] ;
//   Vector4<double>& V4_Kp_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Kp_HeliPsi] ;
//   Vector4<double>& V4_Km_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Km_HeliPsi] ;

//   double thePhiPhiDecayPlaneAngle = decayAngleChi( V4_KsKlKpKm_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi, V4_Ks_HeliPsi, V4_Kl_HeliPsi   );


  Vector4<double>& V4_normKpKmDecHeliKsKlKpKm = theData->FourVecs[enumJpsiGamKsKlKKData::V4_normKpKmDecHeliKsKlKpKm];
  Vector4<double>& V4_normKsKlDecHeliKsKlKpKm = theData->FourVecs[enumJpsiGamKsKlKKData::V4_normKsKlDecHeliKsKlKpKm];

  double cosChi=(V4_normKpKmDecHeliKsKlKpKm.Px()*V4_normKsKlDecHeliKsKlKpKm.Px()
		 +V4_normKpKmDecHeliKsKlKpKm.Py()*V4_normKsKlDecHeliKsKlKpKm.Py()
		 +V4_normKpKmDecHeliKsKlKpKm.Pz()*V4_normKsKlDecHeliKsKlKpKm.Pz())
    / (V4_normKpKmDecHeliKsKlKpKm.P()*V4_normKsKlDecHeliKsKlKpKm.P());

  double chi=acos(fabs(cosChi));
  theChiHisto->Fill(chi*180./TMath::Pi(),weight);   
}

void  JpsiGamKsKlKKHist::fillTuple( TNtuple* theTuple, EvtData* theData, double weight)
{
  
  Vector4<double>& V4_KsKlKpKm_HeliPsi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi]  ;
  Vector4<double>& V4_KsKl_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKl_HeliPsi] ;
  Vector4<double>& V4_KpKm_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliPsi] ;
  Vector4<double>& V4_gamma_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_gamma_HeliPsi] ;
  
  Vector4<double>& V4_KsKl_HeliKsKlKpKm= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKl_HeliKsKlKpKm] ;
  Vector4<double>& V4_KpKm_HeliKsKlKpKm= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KpKm_HeliKsKlKpKm] ;
  Vector4<double>& V4_Ks_HeliKsKl= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Ks_HeliKsKl] ;
  Vector4<double>& V4_Kp_HeliKpKm= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Kp_HeliKpKm] ;
  
  Vector4<double>& V4_Ks_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Ks_HeliPsi] ;
  Vector4<double>& V4_Kl_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Kl_HeliPsi] ;
  Vector4<double>& V4_Kp_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Kp_HeliPsi] ;
  Vector4<double>& V4_Km_HeliPsi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_Km_HeliPsi] ;
  
  
  Vector4<double>& V4_normKpKmDecHeliKsKlKpKm = theData->FourVecs[enumJpsiGamKsKlKKData::V4_normKpKmDecHeliKsKlKpKm];
  Vector4<double>& V4_normKsKlDecHeliKsKlKpKm = theData->FourVecs[enumJpsiGamKsKlKKData::V4_normKsKlDecHeliKsKlKpKm];
  double cosChi=(V4_normKpKmDecHeliKsKlKpKm.Px()*V4_normKsKlDecHeliKsKlKpKm.Px()
		 +V4_normKpKmDecHeliKsKlKpKm.Py()*V4_normKsKlDecHeliKsKlKpKm.Py()
		 +V4_normKpKmDecHeliKsKlKpKm.Pz()*V4_normKsKlDecHeliKsKlKpKm.Pz())
    / (V4_normKpKmDecHeliKsKlKpKm.P()*V4_normKsKlDecHeliKsKlKpKm.P());

  double chi=acos(fabs(cosChi));
  double thePhiPhiDecayPlaneAngle = chi/TMath::Pi()*180.;
  
  //double thePhiPhiDecayPlaneAngle = decayAngleChi( V4_KsKlKpKm_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi, V4_Ks_HeliPsi, V4_Kl_HeliPsi   );
  double testHeli = costDecHeli( V4_KsKlKpKm_HeliPsi+V4_gamma_HeliPsi, V4_Ks_HeliPsi+V4_Kl_HeliPsi+V4_Km_HeliPsi+V4_Kp_HeliPsi, V4_Ks_HeliPsi+V4_Kl_HeliPsi );
  
  theTuple->Fill( V4_KsKlKpKm_HeliPsi.M(),
		  V4_KsKl_HeliPsi.M(),
		  V4_KpKm_HeliPsi.M(),
		  V4_KsKlKpKm_HeliPsi.CosTheta(),
		  V4_gamma_HeliPsi.CosTheta(),
		  V4_KsKl_HeliKsKlKpKm.CosTheta(),
		  V4_KpKm_HeliKsKlKpKm.CosTheta(),
		  V4_Ks_HeliKsKl.CosTheta(),
		  V4_Kp_HeliKpKm.CosTheta(),
		  thePhiPhiDecayPlaneAngle,		  
		  testHeli,
		  weight);
}


double JpsiGamKsKlKKHist::decayAngleChi(const Vector4<double>& v4_p,const Vector4<double>& v4_d1,
		const Vector4<double>& v4_d2,const Vector4<double>& v4_h1,
		const Vector4<double>& v4_h2 ) {
  
  TLorentzVector p4_p(  v4_p.Px(), v4_p.Py(), v4_p.Pz(), v4_p.E() );
  TLorentzVector p4_d1p( v4_d1.Px(), v4_d1.Py(), v4_d1.Pz(), v4_d1.E() );
  TLorentzVector p4_d2p( v4_d2.Px(), v4_d2.Py(), v4_d2.Pz(), v4_d2.E() );
  TLorentzVector p4_h1p( v4_h1.Px(), v4_h1.Py(), v4_h1.Pz(), v4_h1.E() );
  TLorentzVector p4_h2p( v4_h2.Px(), v4_h2.Py(), v4_h2.Pz(), v4_h2.E() );
  
  
  // boost all vectors parent restframe
  
  p4_d1p.Boost( -p4_p.BoostVector() );
  p4_d2p.Boost( -p4_p.BoostVector() );
  p4_h1p.Boost( -p4_p.BoostVector() );
  p4_h2p.Boost( -p4_p.BoostVector() );
  
  
  TVector3 d1_perp,d1_prime,h1_perp;
  TVector3 D;
  
  D=(p4_d1p+p4_d2p).Vect();
  
  d1_perp=p4_d1p.Vect()-(D.Dot(p4_d1p.Vect())/D.Dot(D))*D;
  h1_perp=p4_h1p.Vect()-(D.Dot(p4_h1p.Vect())/D.Dot(D))*D;
  
  // orthogonal to both D and d1_perp
  
  d1_prime=D.Cross(d1_perp);
  
  d1_perp= d1_perp* (1./d1_perp.Mag());
  d1_prime= d1_prime * (1./d1_prime.Mag());
  
  double x,y;
 
 x=d1_perp.Dot(h1_perp);
 y=d1_prime.Dot(h1_perp);
 
 double chi=atan2(y,x);
 
 if (chi<0.0) chi+=2.*TMath::Pi();
 
 return chi;
 
}

