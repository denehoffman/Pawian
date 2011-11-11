#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamHist.hh"
// #include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
// #include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamHist::Psi2SToKpKmPiGamHist(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _cosPsiDataHist(0),
  _cosPsiMcHist(0),
  _cosPsiFittedHist(0),
  _cosChic1Toa0PiDataHist(0),
  _cosChic1Toa0PiMcHist(0),
  _cosChic1Toa0PiFittedHist(0),
  _cosChic1ToK890KDataHist(0),
  _cosChic1ToK890KMcHist(0),
  _cosChic1ToK890KFittedHist(0),
  _cosChic1ToK1400KDataHist(0),
  _cosChic1ToK1400KMcHist(0),
  _cosChic1ToK1400KFittedHist(0),
  _phiChic1ToK1400KDataHist(0),
  _phiChic1ToK1400KMcHist(0),
  _phiChic1ToK1400KFittedHist(0),
  _cosChic1ToKKPiDataHist(0),
  _cosChic1ToKKPiMcHist(0),
  _cosChic1ToKKPiFittedHist(0),
  _phiChic1ToKKPiDataHist(0),
  _phiChic1ToKKPiMcHist(0),
  _phiChic1ToKKPiFittedHist(0),
  _cosK890DataHist(0),
  _cosK890McHist(0),
  _cosK890FittedHist(0),
  _cosK1400DataHist(0),
  _cosK1400McHist(0),
  _cosK1400FittedHist(0),
  _cosa980DataHist(0),
  _cosa980McHist(0),
  _cosa980FittedHist(0),
  _invKpKmDataHist(0),
  _invKpKmMcHist(0),
  _invKpKmFittedHist(0),
  _invKPiDataHist(0),
  _invKPiMcHist(0),
  _invKPiFittedHist(0),
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==theEvtList){
    Alert <<"EvtList* theEvtList is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  initRootStuff();

  const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotCosPsi(_cosPsiDataHist, (*it), 1.);
      plotChic1Toa0Pi(_cosChic1Toa0PiDataHist, (*it), 1.);
      plotChic1ToK890K(_cosChic1ToK890KDataHist, (*it), 1.);
      plotChic1ToK1400K(_cosChic1ToK1400KDataHist,  _phiChic1ToK1400KDataHist, (*it), 1.);
      plotAnglesChic1ToKKPi(_cosChic1ToKKPiDataHist, _phiChic1ToKKPiDataHist, (*it), 1.);
      plotMKpKm(_invKpKmDataHist, (*it), 1.);
      plotMKPi(_invKPiDataHist, (*it), 1.);
      plotCosK890(_cosK890DataHist, (*it), 1.);
      plotCosK1400(_cosK1400DataHist, (*it), 1.);
      plotCosa980(_cosa980DataHist, (*it), 1.);

      //     fillTuple(_dataTuple, (*it), 1.);

      ++it;
    }

  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotCosPsi(_cosPsiMcHist, (*it), 1.);
      plotChic1Toa0Pi(_cosChic1Toa0PiMcHist, (*it), 1.);
      plotChic1ToK890K(_cosChic1ToK890KMcHist, (*it), 1.);
      plotChic1ToK1400K(_cosChic1ToK1400KMcHist,  _phiChic1ToK1400KMcHist, (*it), 1.);
      plotAnglesChic1ToKKPi(_cosChic1ToKKPiMcHist, _phiChic1ToKKPiMcHist, (*it), 1.);
      plotMKpKm(_invKpKmMcHist, (*it), 1.);
      plotMKPi(_invKPiMcHist, (*it), 1.);
      plotCosK890(_cosK890McHist, (*it), 1.);
      plotCosK1400(_cosK1400McHist, (*it), 1.);
      plotCosa980(_cosa980McHist, (*it), 1.);

//      fillTuple(_mcTuple, (*it), 1.);

      ++it;
    }
}

Psi2SToKpKmPiGamHist::Psi2SToKpKmPiGamHist(boost::shared_ptr<AbsLh> thePsi2SToKpKmPiGamLh, fitParams& fitParam) :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _cosPsiDataHist(0),
  _cosPsiMcHist(0),
  _cosPsiFittedHist(0),
  _cosChic1Toa0PiDataHist(0),
  _cosChic1Toa0PiMcHist(0),
  _cosChic1Toa0PiFittedHist(0),
  _cosChic1ToK890KDataHist(0),
  _cosChic1ToK890KMcHist(0),
  _cosChic1ToK890KFittedHist(0),
  _cosChic1ToK1400KDataHist(0),
  _cosChic1ToK1400KMcHist(0),
  _cosChic1ToK1400KFittedHist(0),
  _phiChic1ToK1400KDataHist(0),
  _phiChic1ToK1400KMcHist(0),
  _phiChic1ToK1400KFittedHist(0),
  _cosChic1ToKKPiDataHist(0),
  _cosChic1ToKKPiMcHist(0),
  _cosChic1ToKKPiFittedHist(0),
  _phiChic1ToKKPiDataHist(0),
  _phiChic1ToKKPiMcHist(0),
  _phiChic1ToKKPiFittedHist(0),
  _cosK890DataHist(0),
  _cosK890McHist(0),
  _cosK890FittedHist(0),
  _cosK1400DataHist(0),
  _cosK1400McHist(0),
  _cosK1400FittedHist(0),
  _cosa980DataHist(0),
  _cosa980McHist(0),
  _cosa980FittedHist(0),
  _invKpKmDataHist(0),
  _invKpKmMcHist(0),
  _invKpKmFittedHist(0),
  _invKPiDataHist(0),
  _invKPiMcHist(0),
  _invKPiFittedHist(0),
  _dataTuple(0),
  _mcTuple(0)
{
  if(0==thePsi2SToKpKmPiGamLh){
    Alert <<"Psi2SToKpKmPiGamLh* thePsi2SToKpKmPiGamLh is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  initRootStuff();

  boost::shared_ptr<const EvtDataBaseList> theEvtList=thePsi2SToKpKmPiGamLh->getEventList();
  const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

  std::vector<EvtData*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotCosPsi(_cosPsiDataHist, (*it), 1.);
      plotChic1Toa0Pi(_cosChic1Toa0PiDataHist, (*it), 1.);
      plotChic1ToK890K(_cosChic1ToK890KDataHist, (*it), 1.);
      plotChic1ToK1400K(_cosChic1ToK1400KDataHist,  _phiChic1ToK1400KDataHist, (*it), 1.);
      plotAnglesChic1ToKKPi(_cosChic1ToKKPiDataHist, _phiChic1ToKKPiDataHist, (*it), 1.);
      plotMKpKm(_invKpKmDataHist, (*it), 1.);
      plotMKPi(_invKPiDataHist, (*it), 1.);
      plotCosK890(_cosK890DataHist, (*it), 1.);
      plotCosK1400(_cosK1400DataHist, (*it), 1.);
      plotCosa980(_cosa980DataHist, (*it), 1.); 

      //     fillTuple(_dataTuple, (*it), 1.);

      ++it;
    }

  const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    { 
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotCosPsi(_cosPsiMcHist, (*it), 1.);
      plotChic1Toa0Pi(_cosChic1Toa0PiMcHist, (*it), 1.);
      plotChic1ToK890K(_cosChic1ToK890KMcHist, (*it), 1.);
      plotChic1ToK1400K(_cosChic1ToK1400KMcHist,  _phiChic1ToK1400KMcHist, (*it), 1.);
      plotAnglesChic1ToKKPi(_cosChic1ToKKPiMcHist, _phiChic1ToKKPiMcHist, (*it), 1.);
      plotMKpKm(_invKpKmMcHist, (*it), 1.);
      plotMKPi(_invKPiMcHist, (*it), 1.);
      plotCosK890(_cosK890McHist, (*it), 1.);
      plotCosK1400(_cosK1400McHist, (*it), 1.);
      plotCosa980(_cosa980McHist, (*it), 1.);

      double evtWeight= thePsi2SToKpKmPiGamLh->calcEvtIntensity((*it), fitParam);
      plotDalitz(_dalitzFittedHist, (*it), evtWeight);
      plotCosPsi(_cosPsiFittedHist, (*it), evtWeight);
      plotChic1Toa0Pi(_cosChic1Toa0PiFittedHist, (*it), evtWeight);
      plotChic1ToK890K(_cosChic1ToK890KFittedHist, (*it), evtWeight);
      plotChic1ToK1400K(_cosChic1ToK1400KFittedHist,  _phiChic1ToK1400KFittedHist, (*it), evtWeight);
      plotAnglesChic1ToKKPi(_cosChic1ToKKPiFittedHist, _phiChic1ToKKPiFittedHist, (*it), evtWeight);
      plotCosK890(_cosK890FittedHist, (*it), evtWeight);
      plotCosK1400(_cosK1400FittedHist, (*it), evtWeight);
      plotCosa980(_cosa980FittedHist, (*it), evtWeight);
      plotMKpKm(_invKpKmFittedHist, (*it), evtWeight);
      plotMKPi(_invKPiFittedHist, (*it), evtWeight);

 //      fillTuple(_mcTuple, (*it), evtWeight);
      ++it;
    }

  _cosPsiDataHist->Sumw2();
//   double integralData=_cosPsiDataHist->Integral();
  double integralData=(double) theEvtList->getDataVecs().size();
  Info <<"No of fit data events  " << integralData << endmsg;
 
//   double integralFitted=_cosPsiFittedHist->Integral();


  double integralFitted=(double) theEvtList->getMcVecs().size();
  Info <<"No of fit events " << integralFitted << endmsg; 

  Info <<"scaling factor  " << integralData/integralFitted << endmsg; 

  Info <<"integral fitted hist " << _cosPsiFittedHist->Integral();

  _cosPsiFittedHist->Scale(integralData/integralFitted);
  _dalitzFittedHist->Scale(integralData/integralFitted);
  _cosChic1Toa0PiFittedHist->Scale(integralData/integralFitted);
  _cosChic1ToK890KFittedHist->Scale(integralData/integralFitted);
  _cosChic1ToK1400KFittedHist->Scale(integralData/integralFitted);
  _cosChic1ToKKPiFittedHist->Scale(integralData/integralFitted);
  _phiChic1ToKKPiFittedHist->Scale(integralData/integralFitted);
  _cosK890FittedHist->Scale(integralData/integralFitted);
  _cosK1400FittedHist->Scale(integralData/integralFitted);
  _phiChic1ToK1400KFittedHist->Scale(integralData/integralFitted);
  _cosa980FittedHist->Scale(integralData/integralFitted);
  _invKpKmFittedHist->Scale(integralData/integralFitted);
  _invKPiFittedHist->Scale(integralData/integralFitted);
}

Psi2SToKpKmPiGamHist::~Psi2SToKpKmPiGamHist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void Psi2SToKpKmPiGamHist::initRootStuff()
{
  std::string rootFileName="./Psi2SToKpKmPiGam.root";

  _theTFile=new TFile(rootFileName.c_str(),"recreate");
  _dalitzDataHist= new TH2F("_dalitzDataHist","Dpl K+K- K+#pi^{0} data",90, 0., 12.0, 70, 0., 10.0);
  _dalitzMcHist= new TH2F("_dalitzMcHist","Dpl K+K- K+#pi^{0} MC",90, 0., 12.0, 70, 0., 10.0);
  _dalitzFittedHist= new TH2F("_dalitzFittedHist","Dpl K+K- K+#pi^{0} fit",90, 0., 12.0, 70, 0., 10.0);
  _cosPsiDataHist= new TH1F("_cosPsiDataHist","cos(#Theta) #Psi(2S) data",60, -1., 1.);
  _cosPsiMcHist= new TH1F("_cosPsiMcHist","cos(#Theta) #Psi(2S) MC",60, -1., 1.);
  _cosPsiFittedHist= new TH1F("_cosPsiFittedHist","cos(#Theta) #Psi(2S) fit",60, -1., 1.);
  _cosChic1Toa0PiDataHist= new TH1F("_cosChic1Toa0PiDataHist","cos(#Theta) #Chi_{c1} #rightarrow a_{0}(980) #pi^{0} data",60, -1., 1.);
  _cosChic1Toa0PiMcHist= new TH1F("_cosChic1Toa0PiMcHist","cos(#Theta) #Chi_{c1} #rightarrow a_{0}(980) #pi^{0} Mc",60, -1., 1.);
  _cosChic1Toa0PiFittedHist= new TH1F("_cosChic1Toa0PiFittedHist","cos(#Theta) #Chi_{c1} #rightarrow a_{0}(980) #pi^{0} fit",60, -1., 1.);
  _cosChic1ToK890KDataHist= new TH1F("_cosChic1ToK890KDataHist","cos(#Theta) #chi_{c1} #rightarrow K(890) K data",60, -1., 1.);
  _cosChic1ToK890KMcHist= new TH1F("_cosChic1ToK890KMcHist","cos(#Theta) #chi_{c1} #rightarrow K(890) K Mc",60, -1., 1.);
  _cosChic1ToK890KFittedHist= new TH1F("_cosChic1ToK890KFittedHist","cos(#Theta) #chi_{c1} #rightarrow K(890) K fit",60, -1., 1.);
  _cosChic1ToK1400KDataHist= new TH1F("_cosChic1ToK1400KDataHist","cos(#Theta) #chi_{c1} #rightarrow K(1400) K data",60, -1., 1.);
  _cosChic1ToK1400KMcHist= new TH1F("_cosChic1ToK1400KMcHist","cos(#Theta) #chi_{c1} #rightarrow K(1400) K Mc",60, -1., 1.);
  _cosChic1ToK1400KFittedHist= new TH1F("_cosChic1ToK1400KFittedHist","cos(#Theta) #chi_{c1} #rightarrow K(1400) K fit",60, -1., 1.);

  _cosChic1ToKKPiDataHist= new TH1F("_cosChic1ToKKPiDataHist","cos(#Theta) #chi_{c1} #rightarrow K K #pi data",101, -1., 1.);
  _cosChic1ToKKPiMcHist= new TH1F("_cosChic1ToKKPiMcHist","cos(#Theta) #chi_{c1} #rightarrow K K #pi Mc",101, -1., 1.);
  _cosChic1ToKKPiFittedHist= new TH1F("_cosChic1ToKKPiFittedHist","cos(#Theta) #chi_{c1} #rightarrow K K #pi fit",101, -1., 1.);
  _phiChic1ToKKPiDataHist= new TH1F("_phiChic1ToKKPiDataHist","#Phi #chi_{c1} #rightarrow K K #pi data",101, -TMath::Pi(), TMath::Pi());
  _phiChic1ToKKPiMcHist= new TH1F("_phiChic1ToKKPiMcHist","#Phi #chi_{c1} #rightarrow K K #pi Mc",101, -TMath::Pi(), TMath::Pi());
  _phiChic1ToKKPiFittedHist= new TH1F("_phiChic1ToKKPiFittedHist","#Phi #chi_{c1} #rightarrow K K #pi fit",101, -TMath::Pi(), TMath::Pi());




  _cosK890DataHist= new TH1F("_cosK890DataHist","cos(#Theta) K*(890) data",60, -1., 1.);
  _cosK890McHist= new TH1F("_cosK890McHist","cos(#Theta) K*(890) MC",60, -1., 1.);
  _cosK890FittedHist= new TH1F("_cosK890FittedHist","cos(#Theta) K*(890) fit",60, -1., 1.);

  _cosK1400DataHist= new TH1F("_cosK1400DataHist","cos(#Theta) K*(1400) data",60, -1., 1.);
  _cosK1400McHist= new TH1F("_cosK1400McHist","cos(#Theta) K*(1400) MC",60, -1., 1.);
  _cosK1400FittedHist= new TH1F("_cosK1400FittedHist","cos(#Theta) K*(1400) fit",60, -1., 1.);
  _phiChic1ToK1400KDataHist= new TH1F("_phiChic1ToK1400KDataHist","#phi K*(1400) data",60, -TMath::Pi(), TMath::Pi());
  _phiChic1ToK1400KMcHist= new TH1F("_phiChic1ToK1400KMcHist","#phi K*(1400) MC",60, -TMath::Pi(), TMath::Pi());
  _phiChic1ToK1400KFittedHist= new TH1F("_phiChic1ToK1400KFittedHist","#phi K*(1400) fit",60, -TMath::Pi(), TMath::Pi());


  _cosa980DataHist= new TH1F("_cosa980DataHist","cos(#Theta) a(980) data",60, -1., 1.);
  _cosa980McHist= new TH1F("_cosa980McHist","cos(#Theta) a(980) MC",60, -1., 1.);
  _cosa980FittedHist= new TH1F("_cosa980FittedHist","cos(#Theta) a(980) fit",60, -1., 1.);
  _invKpKmDataHist= new TH1F("_invKpKmDataHist","M_{K+K-} data",120, 0.8, 3.6);
  _invKpKmMcHist= new TH1F("_invKpKmMcHist","M_{K+K-} MC",120, 0.8, 3.6);
  _invKpKmFittedHist= new TH1F("_invKpKmFittedHist","M_{K+K-} fit",120, 0.8, 3.6);
  _invKPiDataHist= new TH1F("_invKPiDataHist","M_{K#pi^{0}} data",140, 0.6, 3.4);
  _invKPiMcHist= new TH1F("_invKPiMcHist","M_{K#pi^{0}} MC",140, 0.6, 3.4);
  _invKPiFittedHist= new TH1F("_invKPiFittedHist","M_{K#pi^{0}} fit",140, 0.6, 3.4);

  _dataTuple=new TNtuple("_dataTuple", "data ntuple", "mKpKm:mKpPi:mKmPi:cosChiGam:cosKstpK:cosKstmK:cosa0Pi:cosKK:cosKpPi:cosKmPi:weight");
  _mcTuple=new TNtuple("_mcTuple", "mc ntuple", "mKpKm:mKpPi:mKmPi:cosChiGam:cosKstpK:cosKstmK:cosa0Pi:cosKK:cosKpPi:cosKmPi:weight");
}

void Psi2SToKpKmPiGamHist::plotDalitz(TH2F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& KpKm_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1];
  Vector4<double>& KpPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  theHisto->Fill(KpKm_HeliChic1_4V.M()*KpKm_HeliChic1_4V.M(), KpPi_HeliChic1_4V.M()*KpPi_HeliChic1_4V.M(),weight);  
}

void Psi2SToKpKmPiGamHist::plotCosPsi(TH1F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& chic1_HeliPsi2S_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::Chic1_HeliPsi2S];
  theHisto->Fill(chic1_HeliPsi2S_4V.CosTheta(), weight);
}

void Psi2SToKpKmPiGamHist::plotChic1Toa0Pi(TH1F* theHisto, EvtData* theData, double weight)
{ 
 Vector4<double>& KpKm_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1];
 if (KpKm_HeliChic1_4V.M()<1.35) theHisto->Fill(KpKm_HeliChic1_4V.CosTheta(), weight); 
}

void Psi2SToKpKmPiGamHist::plotChic1ToK890K(TH1F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& KpPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  if (KpPi_HeliChic1_4V.M()>0.85 && KpPi_HeliChic1_4V.M()<0.95) theHisto->Fill(KpPi_HeliChic1_4V.CosTheta(), weight); 

  Vector4<double>& KmPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];
  if (KmPi_HeliChic1_4V.M()>0.85 && KmPi_HeliChic1_4V.M()<0.95) theHisto->Fill(KmPi_HeliChic1_4V.CosTheta(), weight); 
}

void Psi2SToKpKmPiGamHist::plotChic1ToK1400K(TH1F* costHist, TH1F* phiHist, EvtData* theData, double weight)
{
  Vector4<double> KpPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  if (KpPi_HeliChic1_4V.M()>1.36 && KpPi_HeliChic1_4V.M()<1.48){
    costHist->Fill(KpPi_HeliChic1_4V.CosTheta(), weight);
    phiHist->Fill(KpPi_HeliChic1_4V.Phi(), weight);
  } 

  Vector4<double> KmPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];
  if (KmPi_HeliChic1_4V.M()>1.36 && KmPi_HeliChic1_4V.M()<1.48){
    costHist->Fill(KmPi_HeliChic1_4V.CosTheta(), weight);
    phiHist->Fill(KmPi_HeliChic1_4V.Phi(), weight);
  } 
}

void Psi2SToKpKmPiGamHist::plotAnglesChic1ToKKPi(TH1F* costHist, TH1F* phiHist, EvtData* theData, double weight){
  Vector4<double>& KpKmPi0_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKmPi0_HeliChic1];
  costHist->Fill(KpKmPi0_HeliChic1_4V.CosTheta(), weight);
  phiHist->Fill(KpKmPi0_HeliChic1_4V.Phi(), weight);
}


void Psi2SToKpKmPiGamHist::plotMKpKm(TH1F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& KpKm_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1];
  theHisto->Fill(KpKm_HeliChic1_4V.M(), weight);
}

void Psi2SToKpKmPiGamHist::plotMKPi(TH1F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& KpPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  theHisto->Fill(KpPi_HeliChic1_4V.M(), weight); 
  Vector4<double>& KmPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1]; 
  theHisto->Fill(KmPi_HeliChic1_4V.M(), weight); 
}

void Psi2SToKpKmPiGamHist::plotCosK890(TH1F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& KpPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  Vector4<double>& Kp_HeliKpPi_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::Kp_HeliKpPi];
  if (KpPi_HeliChic1_4V.M()>0.85 && KpPi_HeliChic1_4V.M()<0.95) theHisto->Fill(Kp_HeliKpPi_4V.CosTheta(), weight); 

  Vector4<double>& KmPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];
  Vector4<double>& Km_HeliKmPi_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::Km_HeliKmPi];
  if (KmPi_HeliChic1_4V.M()>0.85 && KmPi_HeliChic1_4V.M()<0.95) theHisto->Fill(Km_HeliKmPi_4V.CosTheta(), weight);
}

void Psi2SToKpKmPiGamHist::plotCosK1400(TH1F* costHist, EvtData* theData, double weight)
{
  Vector4<double>& KpPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  Vector4<double>& Kp_HeliKpPi_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::Kp_HeliKpPi];
  if (KpPi_HeliChic1_4V.M()>1.36 && KpPi_HeliChic1_4V.M()<1.48){
    costHist->Fill(Kp_HeliKpPi_4V.CosTheta(), weight);
//     phiHist->Fill(Kp_HeliKpPi_4V.Phi(), weight);
  } 

  Vector4<double>& KmPi_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];
  Vector4<double>& Km_HeliKmPi_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::Km_HeliKmPi];
  if (KmPi_HeliChic1_4V.M()>1.36 && KmPi_HeliChic1_4V.M()<1.48){
    costHist->Fill(Km_HeliKmPi_4V.CosTheta(), weight);
//     phiHist->Fill(Km_HeliKpPi_4V.Phi(), weight);
  }
}

void Psi2SToKpKmPiGamHist::plotCosa980(TH1F* theHisto, EvtData* theData, double weight)
{
  Vector4<double>& KpKm_HeliChic1_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1];
  Vector4<double>& Kp_HeliKpKm_4V=theData->FourVecs[enumChic1ToKpKmPiGamData::Kp_HeliKpKm];
  if (KpKm_HeliChic1_4V.M()<1.35) theHisto->Fill(Kp_HeliKpKm_4V.CosTheta(), weight); 

}

// void  Psi2SToKpKmPiGamHist::fillTuple( TNtuple* theTuple, const Psi2SToKpKmPiGamEvtData* theData, double weight)
// {
//   Vector4<float> KpKm_HeliChic1_4V=theData->KpKm_HeliChic1_4V;
//   Vector4<float> KpPi_HeliChic1_4V=theData->KpPi_HeliChic1_4V;
//   Vector4<float> KmPi_HeliChic1_4V=theData->KmPi_HeliChic1_4V;
//   Vector4<float> chic1_HeliPsi2S_4V=theData->chic1_HeliPsi2S_4V;
//   Vector4<float> Kp_HeliKpKm_4V=theData->Kp_HeliKpKm_4V;
//   Vector4<float> Kp_HeliKpPi_4V=theData->Kp_HeliKpPi_4V;
//   Vector4<float> Km_HeliKmPi_4V=theData->Km_HeliKmPi_4V;

//   float mKpKm=KpKm_HeliChic1_4V.M();
//   float mKpPi=KpPi_HeliChic1_4V.M();
//   float mKmPi=KmPi_HeliChic1_4V.M();
//   float cosChiGam=chic1_HeliPsi2S_4V.CosTheta();
//   float cosKstpK=KpPi_HeliChic1_4V.CosTheta();
//   float cosKstmK=KmPi_HeliChic1_4V.CosTheta();
//   float cosa0Pi=KpKm_HeliChic1_4V.CosTheta();
//   float cosKK=Kp_HeliKpKm_4V.CosTheta();
//   float cosKpPi=Kp_HeliKpPi_4V.CosTheta();
//   float cosKmPi=Km_HeliKmPi_4V.CosTheta();
//   theTuple->Fill(mKpKm,mKpPi,mKmPi,cosChiGam,cosKstpK,cosKstmK,cosa0Pi, cosKK, cosKpPi, cosKmPi, weight);
// }
