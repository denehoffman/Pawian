#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamHist.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "ErrLogger/ErrLogger.hh"

Psi2STo2K2PiGamHist::Psi2STo2K2PiGamHist(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList) :
  _theTFile(0),
  _cosPsiDataHist(0),
  _cosPsiMcHist(0),
  _cosPsiFittedHist(0),
  _cosK1430DataHist(0),
  _cosK1430McHist(0),
  _cosK1430FittedHist(0),
  _invKKDataHist(0),
  _invKKMcHist(0),
  _invKKFittedHist(0),
  _invKPiDataHist(0),
  _invKPiMcHist(0),
  _invKPiFittedHist(0),
  _invPiPiDataHist(0),
  _invPiPiMcHist(0),
  _invPiPiFittedHist(0),
  _invKPiPiViaK892DataHist(0),
  _invKPiPiViaK892McHist(0),
  _invKPiPiViaK892FittedHist(0),
  _KPivsKPiDataHist2d(0),
  _KPivsKPiMcHist2d(0),
  _KPivsKPiFittedHist2d(0),
  _KKvsPiPiDataHist2d(0),
  _KKvsPiPiMcHist2d(0),
  _KKvsPiPiFittedHist2d(0),
  _KPiPivsPiPiDataHist2d(0),
  _KPiPivsPiPiMcHist2d(0),
  _KPiPivsPiPiFittedHist2d(0),  
  _KKPivsKPiDataHist2d(0),
  _KKPivsKPiMcHist2d(0),
  _KKPivsKPiFittedHist2d(0),
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==theEvtList){
    Alert <<"Psi2SToKpKmPiGamEventList* theEvtList is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  initRootStuff();

  const std::vector<Psi2STo2K2PiGamEvtData*> dataList=theEvtList->getDataVecs();

  std::vector<Psi2STo2K2PiGamEvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotCosPsi(_cosPsiDataHist, (*it), 1.);
      plotInvKK(_invKKDataHist, (*it), 1.);
      plotInvKPi(_invKPiDataHist, (*it), 1.);
      plotInvPiPi(_invPiPiDataHist, (*it), 1.);
      plotInvKPiPi(_invKPiPiViaK892DataHist, (*it), 1., 0.892, 0.03);
      plotKPivsKPi(_KPivsKPiDataHist2d, (*it), 1.);
      plotKKvsPiPi(_KKvsPiPiDataHist2d, (*it), 1.);
      plotCosKst(_cosK892DataHist, (*it), 1., 0.892, 0.03);
      plotCosKst(_cosK1430DataHist, (*it), 1., 1.440, 0.04);
      plotCosKstViaK892(_cosK1430ViaK892DataHist, (*it), 1., 1.430, 0.06);
      plotKPiPivsPiPi(_KPiPivsPiPiDataHist2d, (*it), 1.);
      plotKKPivsKPi(_KKPivsKPiDataHist2d, (*it), 1.);
      writeNTuple(_dataTuple, (*it), 1.);
      ++it;
    }

  const std::vector<Psi2STo2K2PiGamEvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotCosPsi(_cosPsiMcHist, (*it), 1.);
      plotInvKK(_invKKMcHist, (*it), 1.);
      plotInvKPi(_invKPiMcHist, (*it), 1.);
      plotInvPiPi(_invPiPiMcHist, (*it), 1.);
      plotInvKPiPi(_invKPiPiViaK892McHist, (*it), 1., 0.892, 0.03);
      plotKPivsKPi(_KPivsKPiMcHist2d, (*it), 1.);
      plotKKvsPiPi(_KKvsPiPiMcHist2d, (*it), 1.);
      plotCosKst(_cosK892McHist, (*it), 1., 0.892, 0.02);
      plotCosKst(_cosK1430McHist, (*it), 1., 1.440, 0.04);
      plotCosKstViaK892(_cosK1430ViaK892McHist, (*it), 1., 1.430, 0.08);
      plotKPiPivsPiPi(_KPiPivsPiPiMcHist2d, (*it), 1.);
      plotKKPivsKPi(_KKPivsKPiMcHist2d, (*it), 1.);
      writeNTuple(_mcTuple, (*it), 1.);
      ++it;
    }
}

Psi2STo2K2PiGamHist::Psi2STo2K2PiGamHist(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> thePsi2STo2K2PiGamLh, Psi2STo2K2PiGamData::fitParamVal& fitParam) :
  _theTFile(0),
  _cosPsiDataHist(0),
  _cosPsiMcHist(0),
  _cosPsiFittedHist(0),
  _cosK1430DataHist(0),
  _cosK1430McHist(0),
  _cosK1430FittedHist(0),
  _invKKDataHist(0),
  _invKKMcHist(0),
  _invKKFittedHist(0),
  _invKPiDataHist(0),
  _invKPiMcHist(0),
  _invKPiFittedHist(0),
  _invPiPiDataHist(0),
  _invPiPiMcHist(0),
  _invPiPiFittedHist(0),
  _invKPiPiViaK892DataHist(0),
  _invKPiPiViaK892McHist(0),
  _invKPiPiViaK892FittedHist(0),
  _KPivsKPiDataHist2d(0),
  _KPivsKPiMcHist2d(0),
  _KPivsKPiFittedHist2d(0),
  _KKvsPiPiDataHist2d(0),
  _KKvsPiPiMcHist2d(0),
  _KKvsPiPiFittedHist2d(0),
  _KPiPivsPiPiDataHist2d(0),
  _KPiPivsPiPiMcHist2d(0),
  _KPiPivsPiPiFittedHist2d(0),
  _KKPivsKPiDataHist2d(0),
  _KKPivsKPiMcHist2d(0),
  _KKPivsKPiFittedHist2d(0),
  _dataTuple(0),
  _mcTuple(0)
{

  if(0==thePsi2STo2K2PiGamLh){
    Alert <<"Psi2STo2K2PiGamLh* thePsi2STo2K2PiGamLh is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList=thePsi2STo2K2PiGamLh->getEventList();
  const std::vector<Psi2STo2K2PiGamEvtData*> dataList=theEvtList->getDataVecs();

  initRootStuff();

  std::vector<Psi2STo2K2PiGamEvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotCosPsi(_cosPsiDataHist, (*it), 1.);
      plotInvKK(_invKKDataHist, (*it), 1.);
      plotInvKPi(_invKPiDataHist, (*it), 1.);
      plotInvPiPi(_invPiPiDataHist, (*it), 1.);
      plotInvKPiPi(_invKPiPiViaK892DataHist, (*it), 1., 0.892, 0.03);
      plotKPivsKPi(_KPivsKPiDataHist2d, (*it), 1.);
      plotKKvsPiPi(_KKvsPiPiDataHist2d, (*it), 1.);
      plotCosKst(_cosK892DataHist, (*it), 1., 0.892, 0.03);
      plotCosKst(_cosK1430DataHist, (*it), 1., 1.440, 0.04);
      plotCosKstViaK892(_cosK1430ViaK892DataHist, (*it), 1., 1.430, 0.06);
      plotKPiPivsPiPi(_KPiPivsPiPiDataHist2d, (*it), 1.);
      plotKKPivsKPi(_KKPivsKPiDataHist2d, (*it), 1.);
      writeNTuple(_dataTuple, (*it), 1.);
      ++it;
    }

  const std::vector<Psi2STo2K2PiGamEvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotCosPsi(_cosPsiMcHist, (*it), 1.);
      plotInvKK(_invKKMcHist, (*it), 1.);
      plotInvKPi(_invKPiMcHist, (*it), 1.);
      plotInvPiPi(_invPiPiMcHist, (*it), 1.);
      plotInvKPiPi(_invKPiPiViaK892McHist, (*it), 1., 0.892, 0.03);
      plotKPivsKPi(_KPivsKPiMcHist2d, (*it), 1.);
      plotKKvsPiPi(_KKvsPiPiMcHist2d, (*it), 1.);
      plotCosKst(_cosK892McHist, (*it), 1., 0.892, 0.02);
      plotCosKst(_cosK1430McHist, (*it), 1., 1.440, 0.04);
      plotCosKstViaK892(_cosK1430ViaK892McHist, (*it), 1., 1.430, 0.08);
      plotKPiPivsPiPi(_KPiPivsPiPiMcHist2d, (*it), 1.);
      plotKKPivsKPi(_KKPivsKPiMcHist2d, (*it), 1.);

      double evtWeight= thePsi2STo2K2PiGamLh->calcEvtIntensity((*it), fitParam);

      plotCosPsi(_cosPsiFittedHist, (*it), evtWeight);
      plotInvKK(_invKKFittedHist, (*it), evtWeight);
      plotInvKPi(_invKPiFittedHist, (*it), evtWeight);
      plotInvPiPi(_invPiPiFittedHist, (*it), evtWeight);
      plotInvKPiPi(_invKPiPiViaK892FittedHist, (*it), evtWeight, 0.892, 0.03);
      plotKPivsKPi(_KPivsKPiFittedHist2d, (*it), evtWeight);
      plotKKvsPiPi(_KKvsPiPiFittedHist2d, (*it), evtWeight);
      plotCosKst(_cosK892FittedHist, (*it), evtWeight, 0.892, 0.02);
      plotCosKst(_cosK1430FittedHist, (*it), evtWeight, 1.440, 0.04);
      plotCosKstViaK892(_cosK1430ViaK892FittedHist, (*it), evtWeight, 1.430, 0.08);
      plotKPiPivsPiPi(_KPiPivsPiPiFittedHist2d, (*it), evtWeight);
      plotKKPivsKPi(_KKPivsKPiFittedHist2d, (*it), evtWeight);
      writeNTuple(_mcTuple, (*it), evtWeight);
      ++it;
    }

  _cosPsiDataHist->Sumw2();
   double integralData=_cosPsiDataHist->Integral();
//  double integralData=(double) theEvtList->getDataVecs().size();
   Info <<"No of fit data events  " << integralData << endmsg;

   double integralFitted=_cosPsiFittedHist->Integral();
//   double integralFitted=(double) theEvtList->getMcVecs().size();
  Info <<"No of fit events " << integralFitted << endmsg; 

  Info <<"scaling factor  " << integralData/integralFitted << endmsg; 

  Info <<"integral fitted hist " << _cosPsiFittedHist->Integral();

  _cosPsiFittedHist->Scale(integralData/integralFitted);
  _invKKFittedHist->Scale(integralData/integralFitted);
  _invKPiFittedHist->Scale(integralData/integralFitted);
  _invPiPiFittedHist->Scale(integralData/integralFitted);
  _invKPiPiViaK892FittedHist->Scale(integralData/integralFitted);
  _KPivsKPiFittedHist2d->Scale(integralData/integralFitted);
  _KKvsPiPiFittedHist2d->Scale(integralData/integralFitted);
  _cosK892FittedHist->Scale(integralData/integralFitted);
  _cosK1430FittedHist->Scale(integralData/integralFitted);
  _cosK1430ViaK892FittedHist->Scale(integralData/integralFitted);
  _KPiPivsPiPiFittedHist2d->Scale(integralData/integralFitted);
  _KKPivsKPiFittedHist2d->Scale(integralData/integralFitted);
}


Psi2STo2K2PiGamHist::~Psi2STo2K2PiGamHist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void Psi2STo2K2PiGamHist::initRootStuff()
{
  std::string rootFileName="./Psi2STo2K2PiGam.root";

  _theTFile=new TFile(rootFileName.c_str(),"recreate");
  _cosPsiDataHist= new TH1F("_cosPsiDataHist","cos(#Theta) #Psi(2S) data",60, -1., 1.);
  _cosPsiMcHist= new TH1F("_cosPsiMcHist","cos(#Theta) #Psi(2S) MC",60, -1., 1.);
  _cosPsiFittedHist= new TH1F("_cosPsiFittedHist","cos(#Theta) #Psi(2S) fit",60, -1., 1.);
  _cosK892DataHist = new TH1F("_cosK892DataHist","cos(#Theta) K892 dec data",60, -1., 1.);
  _cosK892McHist = new TH1F("_cosK892McHist","cos(#Theta) K892 dec MC",60, -1., 1.);
  _cosK892FittedHist = new TH1F("_cosK892FittedHist","cos(#Theta) K892 dec fit",60, -1., 1.);
  _cosK1430DataHist= new TH1F("_cosK1430DataHist","cos(#Theta) K1430 dec data",60, -1., 1.);
  _cosK1430McHist= new TH1F("_cosK1430McHist","cos(#Theta) K1430 dec MC",60, -1., 1.);
  _cosK1430FittedHist= new TH1F("_cosK1430FittedHist","cos(#Theta) K1430 dec fit",60, -1., 1.);
  _cosK1430ViaK892DataHist= new TH1F("_cosK1430ViaK892DataHist","cos(#Theta) K1430 dec (via K*(892)) data",60, -1., 1.);
  _cosK1430ViaK892McHist= new TH1F("_cosK1430ViaK892McHist","cos(#Theta) K1430 dec (via K*(892)) MC",60, -1., 1.);
  _cosK1430ViaK892FittedHist= new TH1F("_cosK1430ViaK892FittedHist","cos(#Theta) K1430 dec (via K*(892)) fit",60, -1., 1.);
  _invKKDataHist= new TH1F("_invKKDataHist","M_{K+K-} data",120, 0.8, 3.0);
  _invKKMcHist= new TH1F("_invKKMcHist","M_{K+K-} MC",120, 0.8, 3.0);
  _invKKFittedHist= new TH1F("_invKKFittedHist","M_{K+K-} fit",120, 0.8, 3.0);
  _invKPiDataHist= new TH1F("_invKPiDataHist","M_{K+ #pi} data",130, 0.5, 2.8);
  _invKPiMcHist= new TH1F("_invKPiMcHist","M_{K+ #pi} MC",130, 0.5, 2.8);
  _invKPiFittedHist= new TH1F("_invKPiFittedHist","M_{K+ #pi} MC",130, 0.5, 2.8);
  _invPiPiDataHist= new TH1F("_invPiPiDataHist","M_{#pi #pi} data",130, 0.2, 2.5);
  _invPiPiMcHist= new TH1F("_invPiPiMcHist","M_{#pi #pi} MC",130, 0.2, 2.5);
  _invPiPiFittedHist= new TH1F("_invPiPiFittedHist","M_{#pi #pi} fit",130, 0.2, 2.5);
  _invKPiPiViaK892DataHist= new TH1F("_invKPiPiViaK892DataHist","M_{K #pi #pi} (via K*(892)) data",90, 0.9, 3.1);
  _invKPiPiViaK892McHist= new TH1F("_invKPiPiViaK892McHist","M_{K #pi #pi} (via K*(892)) MC",90, 0.9, 3.1);
  _invKPiPiViaK892FittedHist= new TH1F("_invKPiPiViaK892FittedHist","M_{K #pi #pi} (via K*(892)) fit",90, 0.9, 3.1);
  _KPivsKPiDataHist2d= new TH2F("_KPivsKPiDataHist2d","M_{K #pi} vs.M_{K #pi} data", 65, 0.5, 2.8, 65, 0.5, 2.8);
  _KPivsKPiMcHist2d= new TH2F("_KPivsKPiMcHist2d","M_{K #pi} vs.M_{K #pi} MC", 65, 0.5, 2.8, 65, 0.5, 2.8);
  _KPivsKPiFittedHist2d= new TH2F("_KPivsKPiFittedHist2d","M_{K #pi} vs.M_{K #pi} fit", 65, 0.5, 2.8, 65, 0.5, 2.8);
  _KKvsPiPiDataHist2d= new TH2F("_KKvsPiPiDataHist2d","M_{K K} vs.M_{#pi #pi} data", 60, 0.8, 3.0, 60, 0.5, 2.7);
  _KKvsPiPiMcHist2d = new TH2F("_KKvsPiPiMcHist2d","M_{K K} vs.M_{#pi #pi} MC", 60, 0.8, 3.0, 60, 0.5, 2.7);
  _KKvsPiPiFittedHist2d = new TH2F("_KKvsPiPiFittedHist2d","M_{K K} vs.M_{#pi #pi} fit", 60, 0.8, 3.0, 60, 0.5, 2.7);
  _KPiPivsPiPiDataHist2d = new TH2F("_KPiPivsPiPiDataHist2d","M_{K #pi #pi} vs.M_{#pi #pi} data", 60, 0.8, 3.0, 60, 0.3, 2.5);
  _KPiPivsPiPiMcHist2d = new TH2F("_KPiPivsPiPiMcHist2d","M_{K #pi #pi} vs.M_{#pi #pi} MC", 60, 0.8, 3.0, 60, 0.3, 2.5);
  _KPiPivsPiPiFittedHist2d = new TH2F("_KPiPivsPiPiFittedHist2d","M_{K #pi #pi} vs.M_{#pi #pi} fit", 60, 0.8, 3.0, 60, 0.3, 2.5);
  _KKPivsKPiDataHist2d = new TH2F("_KKPivsKPiDataHist2d","M_{K K #pi} vs.M_{K #pi} data", 60, 1.1, 3.3, 60, 0.5, 2.8);
  _KKPivsKPiMcHist2d = new TH2F("_KKPivsKPiMcHist2d","M_{K K #pi} vs.M_{K #pi} MC", 60, 1.1, 3.3, 60, 0.5, 2.8);
  _KKPivsKPiFittedHist2d = new TH2F("_KKPivsKPiFittedHist2d","M_{K K #pi} vs.M_{K #pi} fit", 60, 1.1, 3.3, 60, 0.5, 2.8);

  _dataTuple = new TNtuple("dataTuple","dataTuple","phiKKpipi:costhetaKKpipi:phiK1pi1pi2:costhetaK1pi1pi2:phiK2pi1pi2:costhetaK2pi1pi2:phiK1pi1:costhetaK1pi1:phiK1pi2:costhetaK1pi2:phiK2pi1:costhetaK2pi1:phiK2pi2:costhetaK2pi2:phipi1:costhetapi1:phipi2:costhetapi2:mk1pi1pi2:mk2pi1pi2:mk1pi1:mk1pi2:mk2pi1:mk2pi2:weight");
  _mcTuple = new TNtuple("mcTuple","mcTuple","phiKKpipi:costhetaKKpipi:phiK1pi1pi2:costhetaK1pi1pi2:phiK2pi1pi2:costhetaK2pi1pi2:phiK1pi1:costhetaK1pi1:phiK1pi2:costhetaK1pi2:phiK2pi1:costhetaK2pi1:phiK2pi2:costhetaK2pi2:phipi1:costhetapi1:phipi2:costhetapi2:mk1pi1pi2:mk2pi1pi2:mk1pi1:mk1pi2:mk2pi1:mk2pi2:weight");

}

void Psi2STo2K2PiGamHist::writeNTuple(TNtuple* theTuple, const Psi2STo2K2PiGamEvtData* theData, double weight)
{
  float Phichic0_HeliPsi2S      = theData->chic0_HeliPsi2S_4V.Phi();
  float CosThetachic0_HeliPsi2S = theData->chic0_HeliPsi2S_4V.CosTheta();
  float PhiK1pi1pi2             = theData->KpPiPi_HeliChic0_4V.Phi();
  float CosThetaK1pi1pi2        = theData->KpPiPi_HeliChic0_4V.CosTheta();
  float PhiK2pi1pi2             = theData->KmPiPi_HeliChic0_4V.Phi();
  float CosThetaK2pi1pi2        = theData->KmPiPi_HeliChic0_4V.CosTheta();
  float PhiK1pi1                = theData->KpPi0_HeliKpPi0Pi0_4V.Phi();
  float CosThetaK1pi1           = theData->KpPi0_HeliKpPi0Pi0_4V.CosTheta();
  float PhiK1pi2                = theData->KpPi1_HeliKpPi0Pi0_4V.Phi();
  float CosThetaK1pi2           = theData->KpPi1_HeliKpPi0Pi0_4V.CosTheta();
  float PhiK2pi1                = theData->KmPi0_HeliKmPi0Pi0_4V.Phi();
  float CosThetaK2pi1           = theData->KmPi0_HeliKmPi0Pi0_4V.CosTheta();
  float PhiK2pi2                = theData->KmPi1_HeliKmPi0Pi0_4V.Phi();
  float CosThetaK2pi2           = theData->KmPi1_HeliKmPi0Pi0_4V.CosTheta();
  float Phipi1                  = -theData->Kp_HeliKpPi0_4V.Phi();
  float CosThetapi1             = -theData->Kp_HeliKpPi0_4V.CosTheta();
  float Phipi2                  = -theData->Km_HeliKmPi0_4V.Phi();
  float CosThetapi2             = -theData->Km_HeliKmPi0_4V.CosTheta();
  float mK1pi1pi2               = theData->KpPiPi_HeliChic0_4V.M();
  float mK2pi1pi2               = theData->KmPiPi_HeliChic0_4V.M();
  float mK1pi1                  = theData->KpPi0_HeliChic0_4V.M();
  float mK1pi2                  = theData->KpPi1_HeliChic0_4V.M();
  float MK2pi1                  = theData->KmPi0_HeliChic0_4V.M();
  float MK2pi2                  = theData->KmPi1_HeliChic0_4V.M();
  float evtweight               = weight;

  float fourVectors[50];
  for(int i=0;i<50;i++) {fourVectors[i] = 0;}
  fourVectors[0]  = Phichic0_HeliPsi2S;
  fourVectors[1]  = CosThetachic0_HeliPsi2S;
  fourVectors[2]  = PhiK1pi1pi2;
  fourVectors[3]  = CosThetaK1pi1pi2;
  fourVectors[4]  = PhiK2pi1pi2;
  fourVectors[5]  = CosThetaK2pi1pi2;
  fourVectors[6]  = PhiK1pi1;
  fourVectors[7]  = CosThetaK1pi1;
  fourVectors[8]  = PhiK1pi2;
  fourVectors[9]  = CosThetaK1pi2;
  fourVectors[10] = PhiK2pi1;
  fourVectors[11] = CosThetaK2pi1;
  fourVectors[12] = PhiK2pi2;
  fourVectors[13] = CosThetaK2pi2;
  fourVectors[14] = Phipi1;
  fourVectors[15] = CosThetapi1;
  fourVectors[16] = Phipi2;
  fourVectors[17] = CosThetapi2;
  fourVectors[18] = mK1pi1pi2;
  fourVectors[19] = mK2pi1pi2;
  fourVectors[20] = mK1pi1;
  fourVectors[21] = mK1pi2;
  fourVectors[22] = MK2pi1;
  fourVectors[23] = MK2pi2;
  fourVectors[24] = evtweight;
  //  cout << evtweight << endl;

  theTuple->Fill(fourVectors);
}


void Psi2STo2K2PiGamHist::plotCosPsi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight)
{
  Vector4<float> chic1_HeliPsi2S_4V=theData->chic0_HeliPsi2S_4V;
  theHisto->Fill(chic1_HeliPsi2S_4V.CosTheta(), weight);
}


void Psi2STo2K2PiGamHist::plotCosKst(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width){
  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;
  Vector4<float> Kp_HeliKpPi0_4V=theData->Kp_HeliKpPi0_4V;
  Vector4<float> Kp_HeliKpPi1_4V=theData->Kp_HeliKpPi1_4V;
  Vector4<float> Km_HeliKmPi0_4V=theData->Km_HeliKmPi0_4V;
  Vector4<float> Km_HeliKmPi1_4V=theData->Km_HeliKmPi1_4V;

  if ( fabs(KpPi0_HeliChic0_4V.M()-mass)<width ) theHisto->Fill(Kp_HeliKpPi0_4V.CosTheta());
  if ( fabs(KpPi1_HeliChic0_4V.M()-mass)<width ) theHisto->Fill(Kp_HeliKpPi1_4V.CosTheta());
  if ( fabs(KmPi0_HeliChic0_4V.M()-mass)<width ) theHisto->Fill(Km_HeliKmPi0_4V.CosTheta());
  if ( fabs(KmPi1_HeliChic0_4V.M()-mass)<width ) theHisto->Fill(Km_HeliKmPi1_4V.CosTheta());

}

void Psi2STo2K2PiGamHist::plotCosKstViaK892(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width){
  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  Vector4<float> KpPiPi_HeliChic0_4V=theData->KpPiPi_HeliChic0_4V;
  Vector4<float> KmPiPi_HeliChic0_4V=theData->KmPiPi_HeliChic0_4V;

  Vector4<float> KpPi0_HeliKpPi0Pi0_4V=theData->KpPi0_HeliKpPi0Pi0_4V;
  Vector4<float> KpPi1_HeliKpPi0Pi0_4V=theData->KpPi1_HeliKpPi0Pi0_4V;
  Vector4<float> KmPi0_HeliKmPi0Pi0_4V=theData->KmPi0_HeliKmPi0Pi0_4V;
  Vector4<float> KmPi1_HeliKmPi0Pi0_4V=theData->KmPi1_HeliKmPi0Pi0_4V;

  if (fabs(KpPiPi_HeliChic0_4V.M()-mass)< width){
    if ( fabs(KpPi0_HeliChic0_4V.M()-.892)<0.04 ) theHisto->Fill(KpPi0_HeliKpPi0Pi0_4V.CosTheta());
    if ( fabs(KpPi1_HeliChic0_4V.M()-.892)<0.04 ) theHisto->Fill(KpPi1_HeliKpPi0Pi0_4V.CosTheta());
  }

  if (fabs(KmPiPi_HeliChic0_4V.M()-mass)< width){
    if ( fabs(KmPi0_HeliChic0_4V.M()-.892)<0.04 ) theHisto->Fill(KmPi0_HeliKmPi0Pi0_4V.CosTheta());
    if ( fabs(KmPi1_HeliChic0_4V.M()-.892)<0.04 ) theHisto->Fill(KmPi1_HeliKmPi0Pi0_4V.CosTheta());
  }
}

void Psi2STo2K2PiGamHist::plotInvKK(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){
  Vector4<float> KpKm_HeliChic0_4V=theData->KpKm_HeliChic0_4V;
  theHisto->Fill(KpKm_HeliChic0_4V.M(), weight);
}

void  Psi2STo2K2PiGamHist::plotInvKPi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){
  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  theHisto->Fill(KpPi0_HeliChic0_4V.M(), weight);
  theHisto->Fill(KpPi1_HeliChic0_4V.M(), weight);
  theHisto->Fill(KmPi0_HeliChic0_4V.M(), weight);
  theHisto->Fill(KmPi1_HeliChic0_4V.M(), weight);
}

void  Psi2STo2K2PiGamHist::plotInvPiPi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){
  Vector4<float> PiPi_HeliChic0_4V=theData->PiPi_HeliChic0_4V;
  theHisto->Fill(PiPi_HeliChic0_4V.M(), weight);
}

void Psi2STo2K2PiGamHist::plotInvKPiPi(TH1F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width){
  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  Vector4<float> KpPiPi_HeliChic0_4V=theData->KpPiPi_HeliChic0_4V;
  Vector4<float> KmPiPi_HeliChic0_4V=theData->KmPiPi_HeliChic0_4V;

  if (fabs(KpPi0_HeliChic0_4V.M()-mass)<width || fabs(KpPi1_HeliChic0_4V.M()-mass)<width) theHisto->Fill(KpPiPi_HeliChic0_4V.M(),weight);
  if (fabs(KmPi0_HeliChic0_4V.M()-mass)<width || fabs(KmPi1_HeliChic0_4V.M()-mass)<width) theHisto->Fill(KmPiPi_HeliChic0_4V.M(),weight);
 
}

void Psi2STo2K2PiGamHist::plotKPivsKPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){
  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  theHisto->Fill(KpPi0_HeliChic0_4V.M(), KmPi1_HeliChic0_4V.M(),weight);
  theHisto->Fill(KpPi1_HeliChic0_4V.M(), KmPi0_HeliChic0_4V.M(),weight);
  theHisto->Fill(KmPi0_HeliChic0_4V.M(), KpPi1_HeliChic0_4V.M(),weight);
  theHisto->Fill(KmPi1_HeliChic0_4V.M(), KpPi0_HeliChic0_4V.M(),weight); 
}

void Psi2STo2K2PiGamHist::plotKKvsPiPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){
  Vector4<float> KpKm_HeliChic0_4V=theData->KpKm_HeliChic0_4V;
  Vector4<float> PiPi_HeliChic0_4V=theData->PiPi_HeliChic0_4V;

  theHisto->Fill(KpKm_HeliChic0_4V.M(), PiPi_HeliChic0_4V.M(),weight);
}

void Psi2STo2K2PiGamHist::plotKPiPivsPiPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){
  Vector4<float> KpPiPi_HeliChic0_4V=theData->KpPiPi_HeliChic0_4V;
  Vector4<float> KmPiPi_HeliChic0_4V=theData->KmPiPi_HeliChic0_4V;
  Vector4<float> PiPi_HeliChic0_4V=theData->PiPi_HeliChic0_4V;
  
  theHisto->Fill(KpPiPi_HeliChic0_4V.M(), PiPi_HeliChic0_4V.M(),weight);
  theHisto->Fill(KmPiPi_HeliChic0_4V.M(), PiPi_HeliChic0_4V.M(),weight);
}

void Psi2STo2K2PiGamHist::plotKKPivsKPi(TH2F* theHisto, const Psi2STo2K2PiGamEvtData* theData, double weight){

  Vector4<float> KKPi0_HeliChic0_4V=theData->KKPi0_HeliChic0_4V;
  Vector4<float> KKPi1_HeliChic0_4V=theData->KKPi1_HeliChic0_4V;

  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;  

  theHisto->Fill(KKPi0_HeliChic0_4V.M(), KpPi0_HeliChic0_4V.M(),weight);
  theHisto->Fill(KKPi0_HeliChic0_4V.M(), KmPi0_HeliChic0_4V.M(),weight);
  theHisto->Fill(KKPi1_HeliChic0_4V.M(), KpPi1_HeliChic0_4V.M(),weight);
  theHisto->Fill(KKPi1_HeliChic0_4V.M(), KmPi1_HeliChic0_4V.M(),weight);

}
