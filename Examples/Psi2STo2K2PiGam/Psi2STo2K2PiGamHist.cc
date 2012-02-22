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

Psi2STo2K2PiGamHist::Psi2STo2K2PiGamHist(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::string fileName) :  
  _theTFile(0),
  _dataKey("DataHist"),
  _fittedKey("FittedHist"),
  _mcKey("McHist"),
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==theEvtList){
    Alert <<"Psi2SToKpKmPiGamEventList* theEvtList is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  _keyVec.push_back(_dataKey);
  _keyVec.push_back(_fittedKey);
  _keyVec.push_back(_mcKey);
  initRootStuff(fileName);

  const std::vector<Psi2STo2K2PiGamEvtData*> dataList=theEvtList->getDataVecs();

  std::vector<Psi2STo2K2PiGamEvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      fillHistos(_dataKey, (*it), 1.);
      writeNTuple(_dataTuple, (*it), 1.);
      ++it;
    }

  const std::vector<Psi2STo2K2PiGamEvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    {
      fillHistos(_mcKey, (*it), 1.); 
      writeNTuple(_mcTuple, (*it), 1.);
      ++it;
    }
}

Psi2STo2K2PiGamHist::Psi2STo2K2PiGamHist(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> thePsi2STo2K2PiGamLh, param2K2PiGam& fitParam, const std::string fileName) :
  _theTFile(0),
  _dataKey("DataHist"),
  _fittedKey("FittedHist"),
  _mcKey("McHist"),
  _dataTuple(0),
  _mcTuple(0)
{

  if(0==thePsi2STo2K2PiGamLh){
    Alert <<"Psi2STo2K2PiGamLh* thePsi2STo2K2PiGamLh is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  _keyVec.push_back(_dataKey);
  _keyVec.push_back(_fittedKey);
  _keyVec.push_back(_mcKey);
  
  boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList=thePsi2STo2K2PiGamLh->getEventList();
  const std::vector<Psi2STo2K2PiGamEvtData*> dataList=theEvtList->getDataVecs();

  initRootStuff(fileName);

  std::vector<Psi2STo2K2PiGamEvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    { Vector4<float> chic0_HeliPsi2S_4V=(*it)->chic0_HeliPsi2S_4V;
      double cost_chic0_HeliPsi2S=chic0_HeliPsi2S_4V.CosTheta();
      double theWeight=1.;
//       double theWeight=0.5*(3.*cost_chic0_HeliPsi2S*cost_chic0_HeliPsi2S-1.);
//        double theWeight=5.*(3.*cost_chic0_HeliPsi2S*cost_chic0_HeliPsi2S-1.);
//       double theWeight=sqrt(5.)*(0.5*sqrt(5./2.)*(3.*cost_chic0_HeliPsi2S*cost_chic0_HeliPsi2S-1.)+1);
      fillHistos(_dataKey, (*it), theWeight);
      writeNTuple(_dataTuple, (*it), 1.);
      ++it;
    }

  const std::vector<Psi2STo2K2PiGamEvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    {
      fillHistos(_mcKey, (*it), 1.); 

      double evtWeight= thePsi2STo2K2PiGamLh->calcEvtIntensity((*it), fitParam);
      fillHistos(_fittedKey, (*it), evtWeight);
      writeNTuple(_mcTuple, (*it), evtWeight);
      ++it;
    }

}

Psi2STo2K2PiGamHist::Psi2STo2K2PiGamHist(const param2K2PiGam&, const std::string fileName) :
  _theTFile(0),
  _dataKey("SelectedHist"),
  _fittedKey("FittedHist"),
  _mcKey("UnselectedHist")
{
  _keyVec.push_back(_dataKey);
  _keyVec.push_back(_fittedKey);
  _keyVec.push_back(_mcKey);
  
  initRootStuff(fileName);

}


Psi2STo2K2PiGamHist::~Psi2STo2K2PiGamHist()
{
  //scale the fitted histograms
  std::string cosPsiString="cosPsi"+_dataKey;
  TH1F* cosPsiDataHist=_hist1DMap[cosPsiString];
  double integralData=cosPsiDataHist->Integral();
  Info <<"No of data events  " << integralData << endmsg;

  cosPsiString="cosPsi"+_fittedKey;
  TH1F* cosPsiFittedHist=_hist1DMap[cosPsiString];
  double integralFitted=cosPsiFittedHist->Integral();

  if (integralFitted>0.){ //scale now
    Info <<"No of fit events " << integralFitted << endmsg; 
    
    Info <<"scaling factor  " << integralData/integralFitted << endmsg; 
    
    cosPsiFittedHist->Scale(integralData/integralFitted);
    
    std::string cosK892String="cosK892"+_fittedKey;
    TH1F* cosK892FittedHist=_hist1DMap[cosK892String];
    cosK892FittedHist->Scale(integralData/integralFitted);

    std::string cosK1430Key="cosK1430"+_fittedKey;
    TH1F* cosK1430FittedHist=_hist1DMap[cosK1430Key];
    cosK1430FittedHist->Scale(integralData/integralFitted);

    std::string cosK1430ViaK892Key="cosK1430ViaK892"+_fittedKey;
    TH1F* cosK1430ViaK892FittedHist=_hist1DMap[cosK1430ViaK892Key];
    cosK1430ViaK892FittedHist->Scale(integralData/integralFitted);

    std::string invKKKey="invKK"+_fittedKey;
    TH1F* invKKFittedHist=_hist1DMap[invKKKey];
    invKKFittedHist->Scale(integralData/integralFitted);

    std::string invKPiKey="invKPi"+_fittedKey;
    TH1F* invKPiFittedHist=_hist1DMap[invKPiKey];
    invKPiFittedHist->Scale(integralData/integralFitted);

    std::string invPiPiKey="invPiPi"+_fittedKey;
    TH1F* invPiPiFittedHist=_hist1DMap[invPiPiKey];
    invPiPiFittedHist->Scale(integralData/integralFitted);


    std::string invKPiPiViaK892Key="invKPiPiViaK892"+_fittedKey;
    TH1F* invKPiPiViaK892FittedHist=_hist1DMap[invKPiPiViaK892Key];
    invKPiPiViaK892FittedHist->Scale(integralData/integralFitted);

    std::string KPivsKPiKey="KPivsKPi"+_fittedKey;
    TH2F* KPivsKPiFittedHist=_hist2DMap[KPivsKPiKey];
    KPivsKPiFittedHist->Scale(integralData/integralFitted);

    std::string KKvsPiPiKey="KKvsPiPi"+_fittedKey;
    TH2F* KKvsPiPiFittedHist=_hist2DMap[KKvsPiPiKey];
    KKvsPiPiFittedHist->Scale(integralData/integralFitted);

    std::string KPiPivsPiPiKey="KPiPivsPiPi"+_fittedKey;
    TH2F* KPiPivsPiPiFittedHist=_hist2DMap[KPiPivsPiPiKey];
    KPiPivsPiPiFittedHist->Scale(integralData/integralFitted);

    std::string KKPivsKPiKey="KKPivsKPi"+_fittedKey;
    TH2F* KKPivsKPiFittedHist=_hist2DMap[KKPivsKPiKey];
    KKPivsKPiFittedHist->Scale(integralData/integralFitted);
  }

  _theTFile->Write();
  _theTFile->Close();
}

void Psi2STo2K2PiGamHist::initRootStuff(const std::string& fileName)
{

  _theTFile=new TFile(fileName.c_str(),"recreate");

  initHistMap();

  _dataTuple = new TNtuple("dataTuple","dataTuple","phiKKpipi:costhetaKKpipi:phiK1pi1pi2:costhetaK1pi1pi2:phiK2pi1pi2:costhetaK2pi1pi2:phiK1pi1:costhetaK1pi1:phiK1pi2:costhetaK1pi2:phiK2pi1:costhetaK2pi1:phiK2pi2:costhetaK2pi2:phipi1:costhetapi1:phipi2:costhetapi2:mk1k2pi1:mk1k2pi2:mk1k2:mk1pi1pi2:mk2pi1pi2:mk1pi1:mk1pi2:mk2pi1:mk2pi2:mpipi:costhetaK1ViaK1K2:costhetaKKViaK1K2pi1:costhetaKKViaK1K2pi2:costhetapipiViaK1pipi:costhetapipiViaK2pipi:phipipiViaK1pipi:costhetapiViapipi:phipiViapipi:weight:datatype");

  //  _mcTuple = new TNtuple("mcTuple","mcTuple","phiKKpipi:costhetaKKpipi:phiK1pi1pi2:costhetaK1pi1pi2:phiK2pi1pi2:costhetaK2pi1pi2:phiK1pi1:costhetaK1pi1:phiK1pi2:costhetaK1pi2:phiK2pi1:costhetaK2pi1:phiK2pi2:costhetaK2pi2:phipi1:costhetapi1:phipi2:costhetapi2:mK1K2pi1:mK1K2pi2:mK1K2:mk1pi1pi2:mk2pi1pi2:mk1pi1:mk1pi2:mk2pi1:mk2pi2:mpipi:costhetapipiViaK1pipi:costhetapipiViaK2pipi:phipipiViaK1pipi:costhetapiViapipi:phipiViapipi:weight:datatype");

  _mcTuple = new TNtuple("mcTuple","mcTuple","phiKKpipi:costhetaKKpipi:phiK1pi1pi2:costhetaK1pi1pi2:phiK2pi1pi2:costhetaK2pi1pi2:phiK1pi1:costhetaK1pi1:phiK1pi2:costhetaK1pi2:phiK2pi1:costhetaK2pi1:phiK2pi2:costhetaK2pi2:phipi1:costhetapi1:phipi2:costhetapi2:mk1k2pi1:mk1k2pi2:mk1k2:mk1pi1pi2:mk2pi1pi2:mk1pi1:mk1pi2:mk2pi1:mk2pi2:mpipi:costhetaK1ViaK1K2:costhetaKKFromK1K2pi1:costhetaKKFromK1K2pi2:costhetapipiViaK1pipi:costhetapipiViaK2pipi:phipipiViaK1pipi:costhetapiViapipi:phipiViapipi:weight:datatype");
}

void Psi2STo2K2PiGamHist::initHistMap(){
  std::vector<std::string>::const_iterator itKey;
  TH1F* tmpHist;
  for (itKey=_keyVec.begin(); itKey!=_keyVec.end(); ++itKey){
    std::string cosPsiString="cosPsi"+(*itKey);
    std::string cosPsiComment="cos(#Theta) #Psi(2S) "+(*itKey);
    tmpHist= new TH1F(cosPsiString.c_str(),cosPsiComment.c_str(),60, -1., 1.);
    _hist1DMap[cosPsiString]= tmpHist;

    std::string cosK892String="cosK892"+(*itKey);
    std::string cosK892Comment="cos(#Theta) K892 dec "+(*itKey);
    tmpHist = new TH1F(cosK892String.c_str(),cosK892Comment.c_str(),60, -1., 1.);
    _hist1DMap[cosK892String] = tmpHist;

    std::string cosK1430String="cosK1430"+(*itKey);
    std::string cosK1430Comment="cos(#Theta) K1430 dec "+(*itKey);
    tmpHist = new TH1F(cosK1430String.c_str(), cosK1430Comment.c_str() ,60, -1., 1.);
    _hist1DMap[cosK1430String] = tmpHist;

    std::string cosK1430ViaK892String="cosK1430ViaK892"+(*itKey);
    std::string cosK1430ViaK892Comment="cos(#Theta) K1430 dec (via K*(892)) "+(*itKey);
    tmpHist = new TH1F(cosK1430ViaK892String.c_str(), cosK1430ViaK892Comment.c_str() ,60, -1., 1.);
    _hist1DMap[cosK1430ViaK892String] = tmpHist;

    std::string invKKString="invKK"+(*itKey);
    std::string invKKComment="M_{K+K-} "+(*itKey);
    tmpHist = new TH1F(invKKString.c_str(), invKKComment.c_str() ,120, 0.9, 3.3);
    _hist1DMap[invKKString] = tmpHist;

    std::string invKPiString="invKPi"+(*itKey);
    std::string invKPiComment="M_{K+ #pi} "+(*itKey);
    tmpHist = new TH1F(invKPiString.c_str(), invKPiComment.c_str() ,120, 0.5, 2.9);
    _hist1DMap[invKPiString] = tmpHist;

    std::string invPiPiString="invPiPi"+(*itKey);
    std::string invPiPiComment="M_{#pi #pi} "+(*itKey);
    tmpHist = new TH1F(invPiPiString.c_str(), invPiPiComment.c_str() ,120, 0.2, 2.5);
    _hist1DMap[invPiPiString] = tmpHist;

    std::string invKPiPiViaK892String="invKPiPiViaK892"+(*itKey);
    std::string invKPiPiViaK892Comment="M_{K #pi #pi} (via K*892) "+(*itKey);
    tmpHist = new TH1F(invKPiPiViaK892String.c_str(), invKPiPiViaK892Comment.c_str() ,120, 0.8, 3.2);
    _hist1DMap[invKPiPiViaK892String] = tmpHist;

    std::string KPivsKPiString="KPivsKPi"+(*itKey);
    std::string KPivsKPiComment="M_{K #pi} vs M_{K #pi} "+(*itKey);
    TH2F* tmpHist2D=new TH2F(KPivsKPiString.c_str(), KPivsKPiComment.c_str(), 48, 0.5, 2.9, 48, 0.5, 2.9);
    _hist2DMap[KPivsKPiString] = tmpHist2D;

    std::string KKvsPiPiString="KKvsPiPi"+(*itKey);
    std::string KKvsPiPiComment="M_{K K} vs M_{#pi #pi} "+(*itKey);
    tmpHist2D=new TH2F(KKvsPiPiString.c_str(), KKvsPiPiComment.c_str(), 48, 0.9, 3.3, 46, 0.2, 2.5);
    _hist2DMap[KKvsPiPiString] = tmpHist2D;

    std::string KPiPivsPiPiString="KPiPivsPiPi"+(*itKey);
    std::string KPiPivsPiPiComment="M_{K #pi #pi} vs.M_{#pi #pi} "+(*itKey);
    tmpHist2D=new TH2F(KPiPivsPiPiString.c_str(), KPiPivsPiPiComment.c_str(), 46, 0.2, 2.5, 48, 0.8, 3.2);
    _hist2DMap[KPiPivsPiPiString] = tmpHist2D;

    std::string KKPivsKPiString="KKPivsKPi"+(*itKey);
    std::string KKPivsKPiComment="M_{K K #pi} vs.M_{K #pi} "+(*itKey);
    tmpHist2D=new TH2F(KKPivsKPiString.c_str(), KKPivsKPiComment.c_str(), 48, 0.5, 2.9, 48, 1.1, 3.5);
    _hist2DMap[KKPivsKPiString] = tmpHist2D;


    if ((*itKey) ==_dataKey && _dataKey=="DataHist"){

      std::map <std::string, TH1F* >::const_iterator it1d;
      for (it1d=_hist1DMap.begin(); it1d!=_hist1DMap.end(); ++it1d){
	it1d->second->Sumw2();
      }

    }

  }
}

void Psi2STo2K2PiGamHist::fillHistos(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){
      plotCosPsi(theKey, theData, weight);

      std::string costK892Key="cosK892"+theKey;
      plotCosKst(costK892Key, theData, weight, 0.892, 0.02);

      std::string cosK1430Key="cosK1430"+theKey;
      plotCosKst(cosK1430Key, theData, weight, 1.440, 0.04);

      plotCosKstViaK892(theKey, theData, weight, 1.430, 0.08);
      plotInvKK(theKey, theData, weight);
      plotInvKPi(theKey, theData, weight);
      plotInvPiPi(theKey, theData, weight);
      plotInvKPiPi(theKey, theData, weight, 0.892, 0.03);

      plotKPivsKPi(theKey, theData, weight);
      plotKKvsPiPi(theKey, theData, weight);
      plotKPiPivsPiPi(theKey, theData, weight);
      plotKKPivsKPi(theKey, theData, weight);
}

void Psi2STo2K2PiGamHist::writeNTuple(TNtuple* theTuple, const Psi2STo2K2PiGamEvtData* theData, double weight)
{
  float datatype = 1.;  // 1 = real data, 2 = pwamc

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
  float Phipi2                  = -theData->Kp_HeliKpPi1_4V.Phi();
  float CosThetapi2             = -theData->Kp_HeliKpPi1_4V.CosTheta();
  float mK1K2pi1                = theData->KKPi0_HeliChic0_4V.M(); 
  float mK1K2pi2                = theData->KKPi1_HeliChic0_4V.M();  
  float mK1K2                   = theData->KpKm_HeliChic0_4V.M(); 
  float mK1pi1pi2               = theData->KpPiPi_HeliChic0_4V.M();
  float mK2pi1pi2               = theData->KmPiPi_HeliChic0_4V.M();
  float mK1pi1                  = theData->KpPi0_HeliChic0_4V.M();
  float mK1pi2                  = theData->KpPi1_HeliChic0_4V.M();
  float mK2pi1                  = theData->KmPi0_HeliChic0_4V.M();
  float mK2pi2                  = theData->KmPi1_HeliChic0_4V.M();

  float mpipi                   = theData->PiPi_HeliChic0_4V.M();
  float CosThetaK1FromK1K2      = theData->Km_HeliKmKp_4V.CosTheta();
  float CosThetaKKFromK1K2Pi1   = theData->KK_HeliKKPi0_4V.CosTheta(); 
  float CosThetaKKFromK1K2Pi2   = theData->KK_HeliKKPi1_4V.CosTheta(); 
  float CosThetaPiPiFromK1PiPi  = theData->PiPi_HeliKpPi0Pi0_4V.CosTheta();
  float CosThetaPiPiFromK2PiPi  = theData->PiPi_HeliKmPi0Pi0_4V.CosTheta();
  float PhiPiPiFromK1PiPi  = theData->PiPi_HeliKpPi0Pi0_4V.Phi();

  float CosThetaPiFromPiPi      = theData->Pi0_HeliPi0Pi0_4V.CosTheta();
  float PhiPiFromPiPi      = theData->Pi0_HeliPi0Pi0_4V.Phi();

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
  fourVectors[18] = mK1K2pi1;
  fourVectors[19] = mK1K2pi2;
  fourVectors[20] = mK1K2;
  fourVectors[21] = mK1pi1pi2;
  fourVectors[22] = mK2pi1pi2;
  fourVectors[23] = mK1pi1;
  fourVectors[24] = mK1pi2;
  fourVectors[25] = mK2pi1;
  fourVectors[26] = mK2pi2;
  fourVectors[27] = mpipi;
  fourVectors[28] = CosThetaK1FromK1K2;
  fourVectors[29] = CosThetaKKFromK1K2Pi1;
  fourVectors[30] = CosThetaKKFromK1K2Pi2;
  fourVectors[31] = CosThetaPiPiFromK1PiPi;
  fourVectors[32] = CosThetaPiPiFromK2PiPi;
  fourVectors[33] = PhiPiPiFromK1PiPi;
  fourVectors[34] = CosThetaPiFromPiPi;
  fourVectors[35] = PhiPiFromPiPi;

  fourVectors[36] = evtweight;
  fourVectors[37] = datatype;

  //  cout << evtweight << endl;

  theTuple->Fill(fourVectors);
}


void Psi2STo2K2PiGamHist::plotCosPsi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){ 
  std::string cosPsiString="cosPsi"+theKey;
  TH1F* theHisto=_hist1DMap[cosPsiString];
  Vector4<float> chic1_HeliPsi2S_4V=theData->chic0_HeliPsi2S_4V;
  theHisto->Fill(chic1_HeliPsi2S_4V.CosTheta(), weight);
}

void Psi2STo2K2PiGamHist::plotCosKst(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width){
  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;
  Vector4<float> Kp_HeliKpPi0_4V=theData->Kp_HeliKpPi0_4V;
  Vector4<float> Kp_HeliKpPi1_4V=theData->Kp_HeliKpPi1_4V;
  Vector4<float> Km_HeliKmPi0_4V=theData->Km_HeliKmPi0_4V;
  Vector4<float> Km_HeliKmPi1_4V=theData->Km_HeliKmPi1_4V;
  TH1F* theHisto=_hist1DMap[theKey];
  if ( fabs(KpPi0_HeliChic0_4V.M()-mass) < width ) theHisto->Fill(Kp_HeliKpPi0_4V.CosTheta());
  if ( fabs(KpPi1_HeliChic0_4V.M()-mass) < width ) theHisto->Fill(Kp_HeliKpPi1_4V.CosTheta());
  if ( fabs(KmPi0_HeliChic0_4V.M()-mass) < width ) theHisto->Fill(Km_HeliKmPi0_4V.CosTheta());
  if ( fabs(KmPi1_HeliChic0_4V.M()-mass) < width ) theHisto->Fill(Km_HeliKmPi1_4V.CosTheta());

}


void Psi2STo2K2PiGamHist::plotCosKstViaK892(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width){ 

  std::string cosPsiString="cosK1430ViaK892"+theKey;
  TH1F* theHisto=_hist1DMap[cosPsiString];  

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

void Psi2STo2K2PiGamHist::plotInvKK(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){

  std::string invKKKString="invKK"+theKey;
  TH1F* theHisto=_hist1DMap[invKKKString];
  Vector4<float> KpKm_HeliChic0_4V=theData->KpKm_HeliChic0_4V;
  theHisto->Fill(KpKm_HeliChic0_4V.M(), weight); 
}

void  Psi2STo2K2PiGamHist::plotInvKPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){
  std::string invKPiString="invKPi"+theKey;
  TH1F* theHisto=_hist1DMap[invKPiString];

  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  theHisto->Fill(KpPi0_HeliChic0_4V.M(), weight);
  theHisto->Fill(KpPi1_HeliChic0_4V.M(), weight);
  theHisto->Fill(KmPi0_HeliChic0_4V.M(), weight);
  theHisto->Fill(KmPi1_HeliChic0_4V.M(), weight);
}



void  Psi2STo2K2PiGamHist::plotInvPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){
  std::string invPiPiString="invPiPi"+theKey;
  TH1F* theHisto=_hist1DMap[invPiPiString];
  Vector4<float> PiPi_HeliChic0_4V=theData->PiPi_HeliChic0_4V;
  theHisto->Fill(PiPi_HeliChic0_4V.M(), weight);
}



void Psi2STo2K2PiGamHist::plotInvKPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight, double mass, double width){

  std::string invKPiPiViaK892String="invKPiPiViaK892"+theKey;
  TH1F* theHisto=_hist1DMap[invKPiPiViaK892String];

  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  Vector4<float> KpPiPi_HeliChic0_4V=theData->KpPiPi_HeliChic0_4V;
  Vector4<float> KmPiPi_HeliChic0_4V=theData->KmPiPi_HeliChic0_4V;

  if (fabs(KpPi0_HeliChic0_4V.M()-mass)<width || fabs(KpPi1_HeliChic0_4V.M()-mass)<width) theHisto->Fill(KpPiPi_HeliChic0_4V.M(),weight);
  if (fabs(KmPi0_HeliChic0_4V.M()-mass)<width || fabs(KmPi1_HeliChic0_4V.M()-mass)<width) theHisto->Fill(KmPiPi_HeliChic0_4V.M(),weight);
 
}


void Psi2STo2K2PiGamHist::plotKPivsKPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){

  std::string KPivsKPiString="KPivsKPi"+theKey;
  TH2F* theHisto=_hist2DMap[KPivsKPiString];

  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;

  theHisto->Fill(KpPi0_HeliChic0_4V.M(), KmPi1_HeliChic0_4V.M(),weight);
  theHisto->Fill(KpPi1_HeliChic0_4V.M(), KmPi0_HeliChic0_4V.M(),weight);
  theHisto->Fill(KmPi0_HeliChic0_4V.M(), KpPi1_HeliChic0_4V.M(),weight);
  theHisto->Fill(KmPi1_HeliChic0_4V.M(), KpPi0_HeliChic0_4V.M(),weight); 
}



void Psi2STo2K2PiGamHist::plotKKvsPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){

  std::string KKvsPiPiString="KKvsPiPi"+theKey;
  TH2F* theHisto=_hist2DMap[KKvsPiPiString];

  Vector4<float> KpKm_HeliChic0_4V=theData->KpKm_HeliChic0_4V;
  Vector4<float> PiPi_HeliChic0_4V=theData->PiPi_HeliChic0_4V;

  theHisto->Fill(KpKm_HeliChic0_4V.M(), PiPi_HeliChic0_4V.M(),weight);
}


void Psi2STo2K2PiGamHist::plotKPiPivsPiPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){

  std::string KPiPivsPiPiString="KPiPivsPiPi"+theKey;
  TH2F* theHisto=_hist2DMap[KPiPivsPiPiString];
  
  Vector4<float> KpPiPi_HeliChic0_4V=theData->KpPiPi_HeliChic0_4V;
  Vector4<float> KmPiPi_HeliChic0_4V=theData->KmPiPi_HeliChic0_4V;
  Vector4<float> PiPi_HeliChic0_4V=theData->PiPi_HeliChic0_4V;

  theHisto->Fill(PiPi_HeliChic0_4V.M(), KpPiPi_HeliChic0_4V.M(), weight);
  theHisto->Fill(PiPi_HeliChic0_4V.M(), KmPiPi_HeliChic0_4V.M(), weight);
}




void Psi2STo2K2PiGamHist::plotKKPivsKPi(const std::string& theKey, const Psi2STo2K2PiGamEvtData* theData, double weight){

  std::string KKPivsKPiString="KKPivsKPi"+theKey;
  TH2F* theHisto=_hist2DMap[KKPivsKPiString];

  Vector4<float> KKPi0_HeliChic0_4V=theData->KKPi0_HeliChic0_4V;
  Vector4<float> KKPi1_HeliChic0_4V=theData->KKPi1_HeliChic0_4V;

  Vector4<float> KpPi0_HeliChic0_4V=theData->KpPi0_HeliChic0_4V;
  Vector4<float> KpPi1_HeliChic0_4V=theData->KpPi1_HeliChic0_4V;
  Vector4<float> KmPi0_HeliChic0_4V=theData->KmPi0_HeliChic0_4V;
  Vector4<float> KmPi1_HeliChic0_4V=theData->KmPi1_HeliChic0_4V;  

  theHisto->Fill(KpPi0_HeliChic0_4V.M(), KKPi0_HeliChic0_4V.M(), weight);
  theHisto->Fill(KmPi0_HeliChic0_4V.M(), KKPi0_HeliChic0_4V.M(), weight);
  theHisto->Fill(KpPi1_HeliChic0_4V.M(), KKPi1_HeliChic0_4V.M(), weight);
  theHisto->Fill(KmPi1_HeliChic0_4V.M(), KKPi1_HeliChic0_4V.M(), weight);

}


