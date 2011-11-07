#include <getopt.h>
#include <fstream>
#include <sstream>
#include "Examples/pbarpToOmegaPiLS/OmegaPiHistLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "ErrLogger/ErrLogger.hh"

OmegaPiHistLS::OmegaPiHistLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, const std::string &thePathToRootFile) :
  _theTFile(0),
  _cosOmegaHeliDataHist(0),
  _cosOmegaHeliMcHist(0),
  _cosOmegaHeliFittedHist(0),
  _cosOmegaAccCorHist(0),
  _cosOmegaFittedAccCorHist(0),
  _cosPi0FromOmegaDataHeli(0),
  _cosPi0FromOmegaMcHeli(0),
  _cosPi0FromOmegaFittedHeli(0),
  _cosPi0FromOmegaAccCorHeli(0),
  _cosPi0FromOmegaFittedAccCorHeli(0),
  _treimanYangDataHist(0),
  _treimanYangMcHist(0),
  _treimanYangFittedHist(0),
  _treimanYangAccCorHist(0),
  _cosPi0FromOmegaDataHeli1(0),
  _thetaPhiPi0FromOmegaDataHeli(0),
  _thetaPhiPi0FromOmegaMcHeli(0),
  _thetaPhiPi0FromOmegaFittedHeli(0),
  _thetaPhiPi0FromOmegaAccCorHeli(0),
  _thetaPhiPi0FromOmegaFittedAccCorHeli(0),
  _prodDecThetaPhiPi0FromOmegaDataHeli(0),
  _prodDecThetaPhiPi0FromOmegaMcHeli(0),
  _prodDecThetaPhiPi0FromOmegaFittedHeli(0)
{
  if(0==theEvtList){
    Alert <<"OmegaPiEventList* theEvtList is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  _lmax=theEvtList->lMax();
  _pbarmom=theEvtList->pbarMom();

  initRootStuff(thePathToRootFile);

  _cosOmegaHeliDataHist->Sumw2();
  _cosOmegaHeliMcHist->Sumw2();
  _cosPi0FromOmegaDataHeli->Sumw2();
  _cosPi0FromOmegaMcHeli->Sumw2();
  _treimanYangDataHist->Sumw2();
  _thetaPhiPi0FromOmegaDataHeli->Sumw2();
  _thetaPhiPi0FromOmegaMcHeli->Sumw2();
  _prodDecThetaPhiPi0FromOmegaDataHeli->Sumw2();
  _prodDecThetaPhiPi0FromOmegaMcHeli->Sumw2();

  const std::vector<OmPiEvtDataLS*> dataList=theEvtList->getDataVecs();
  std::vector<OmPiEvtDataLS*>::const_iterator it=dataList.begin();
  while(it!=dataList.end())
  {
      double dataEventWeight = ((*it)->eventWeight);
      plotCosOmegaHeli(_cosOmegaHeliDataHist, (*it), dataEventWeight);
      plotCosPi0FromOmegaHeli(_cosPi0FromOmegaDataHeli, (*it), dataEventWeight);
      plotTreimanYang(_treimanYangDataHist, (*it), dataEventWeight);
      plotCosPi0FromOmegaHeli1(_cosPi0FromOmegaDataHeli1, (*it), 1.);
      plotThetaPhiPi0FromOmegaHeli(_thetaPhiPi0FromOmegaDataHeli,(*it), dataEventWeight);                                
      plotProdDecThetaPhiPi0FromOmegaHeli(_prodDecThetaPhiPi0FromOmegaDataHeli, (*it), dataEventWeight); 
      ++it;
  }

  const std::vector<OmPiEvtDataLS*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
  {
      plotCosOmegaHeli(_cosOmegaHeliMcHist, (*it), 1.);
      plotCosPi0FromOmegaHeli(_cosPi0FromOmegaMcHeli, (*it), 1.);
      plotTreimanYang(_treimanYangMcHist, (*it), 1.);
      plotThetaPhiPi0FromOmegaHeli(_thetaPhiPi0FromOmegaMcHeli,(*it),1.);
      plotProdDecThetaPhiPi0FromOmegaHeli(_prodDecThetaPhiPi0FromOmegaMcHeli, (*it), 1.);
      ++it;
  }

  _cosOmegaAccCorHist=doAccCor(_cosOmegaHeliDataHist, _cosOmegaHeliMcHist, "_cosOmegaAccCorHist", "cos(#Theta) #omega heli acc cor");
  _cosPi0FromOmegaAccCorHeli=doAccCor(_cosPi0FromOmegaDataHeli, _cosPi0FromOmegaMcHeli, "_cosPi0FromOmegaAccCorHeli", "cos(#Theta) #pi^0 heli (#omega decay) acc cor");
  _treimanYangAccCorHist=doAccCor(_treimanYangDataHist, _treimanYangMcHist, "_treimanYangAccCorHist", "_treimanYangAccCorHist");
  _thetaPhiPi0FromOmegaAccCorHeli=doAccCor(_thetaPhiPi0FromOmegaDataHeli, _thetaPhiPi0FromOmegaMcHeli,"_thetaPhiPi0FromOmegaAccCorHeli", "_thetaPhiPi0FromOmegaAccCorHeli");
 
}

OmegaPiHistLS::OmegaPiHistLS(boost::shared_ptr<AbsOmegaPiLhLS> omegaPiLh, OmegaPiDataLS::fitParamVal& fitParam, const std::string &thePathToRootFile) :
  _theTFile(0),
  _cosOmegaHeliDataHist(0),
  _cosOmegaHeliMcHist(0),
  _cosOmegaHeliFittedHist(0),
  _cosOmegaAccCorHist(0),
  _cosOmegaFittedAccCorHist(0),
  _cosPi0FromOmegaDataHeli(0),
  _cosPi0FromOmegaMcHeli(0),
  _cosPi0FromOmegaFittedHeli(0),
  _cosPi0FromOmegaAccCorHeli(0),
  _cosPi0FromOmegaFittedAccCorHeli(0),
  _treimanYangDataHist(0),
  _treimanYangMcHist(0),
  _treimanYangFittedHist(0),
  _treimanYangAccCorHist(0),
  _cosPi0FromOmegaDataHeli1(0),
  _thetaPhiPi0FromOmegaDataHeli(0),
  _thetaPhiPi0FromOmegaMcHeli(0),
  _thetaPhiPi0FromOmegaFittedHeli(0),
  _thetaPhiPi0FromOmegaAccCorHeli(0),
  _thetaPhiPi0FromOmegaFittedAccCorHeli(0),
  _prodDecThetaPhiPi0FromOmegaDataHeli(0),
  _prodDecThetaPhiPi0FromOmegaMcHeli(0),
  _prodDecThetaPhiPi0FromOmegaFittedHeli(0)
{
  if(0==omegaPiLh){
    Alert <<"OmegaPiLh* omegaPiLh is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList = omegaPiLh->getEventList();  
  if(0==theEvtList){
    Alert <<"OmegaPiEventList* theEvtList is a 0 pointer !!!!" << endmsg;
    exit(1);
  }

  _lmax=theEvtList->lMax();
  _pbarmom=theEvtList->pbarMom();

  initRootStuff(thePathToRootFile);

  _cosOmegaHeliDataHist->Sumw2();
  _cosOmegaHeliMcHist->Sumw2();
  _cosPi0FromOmegaDataHeli->Sumw2();
  _cosPi0FromOmegaMcHeli->Sumw2();
  _treimanYangDataHist->Sumw2();
  _thetaPhiPi0FromOmegaDataHeli->Sumw2();
  _thetaPhiPi0FromOmegaMcHeli->Sumw2();
  _prodDecThetaPhiPi0FromOmegaDataHeli->Sumw2();
  _prodDecThetaPhiPi0FromOmegaMcHeli->Sumw2();

  std::vector<OmPiEvtDataLS*> dataList=theEvtList->getDataVecs();

  std::vector<OmPiEvtDataLS*>::iterator it=dataList.begin();
  while(it!=dataList.end())
    {
      double dataEventWeight = ((*it)->eventWeight);
      plotCosOmegaHeli(_cosOmegaHeliDataHist, (*it), dataEventWeight);
      plotCosPi0FromOmegaHeli(_cosPi0FromOmegaDataHeli, (*it), dataEventWeight);
      plotTreimanYang(_treimanYangDataHist, (*it), dataEventWeight);
      plotCosPi0FromOmegaHeli1(_cosPi0FromOmegaDataHeli1, (*it), 1.);
      plotThetaPhiPi0FromOmegaHeli(_thetaPhiPi0FromOmegaDataHeli,(*it), dataEventWeight);                                
      plotProdDecThetaPhiPi0FromOmegaHeli(_prodDecThetaPhiPi0FromOmegaDataHeli, (*it), dataEventWeight); 

      ++it;
    }

  std::vector<OmPiEvtDataLS*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    {
      plotCosOmegaHeli(_cosOmegaHeliMcHist, (*it), 1.);
      plotCosPi0FromOmegaHeli(_cosPi0FromOmegaMcHeli, (*it), 1.);
      plotTreimanYang(_treimanYangMcHist, (*it), 1.);
      plotThetaPhiPi0FromOmegaHeli(_thetaPhiPi0FromOmegaMcHeli,(*it),1.);
      plotProdDecThetaPhiPi0FromOmegaHeli(_prodDecThetaPhiPi0FromOmegaMcHeli, (*it), 1.);
      double evtWeight= omegaPiLh->calcEvtIntensity((*it), fitParam);
      plotCosOmegaHeli(_cosOmegaHeliFittedHist, (*it), evtWeight);
      plotCosPi0FromOmegaHeli(_cosPi0FromOmegaFittedHeli, (*it), evtWeight);
      plotTreimanYang(_treimanYangFittedHist, (*it), evtWeight);
      plotThetaPhiPi0FromOmegaHeli(_thetaPhiPi0FromOmegaFittedHeli, (*it), evtWeight);
      plotProdDecThetaPhiPi0FromOmegaHeli(_prodDecThetaPhiPi0FromOmegaFittedHeli, (*it), evtWeight);
      ++it;
    }

  _cosOmegaAccCorHist=doAccCor(_cosOmegaHeliDataHist, _cosOmegaHeliMcHist, "_cosOmegaAccCorHist", "cos(#Theta) #omega heli acc cor");
  _cosPi0FromOmegaAccCorHeli=doAccCor(_cosPi0FromOmegaDataHeli, _cosPi0FromOmegaMcHeli, "_cosPi0FromOmegaAccCorHeli", "cos(#Theta) #pi^0 heli (#omega decay) acc cor");
  _treimanYangAccCorHist=doAccCor(_treimanYangDataHist, _treimanYangMcHist, "_treimanYangAccCorHist", "_treimanYangAccCorHist");
  _thetaPhiPi0FromOmegaAccCorHeli=doAccCor(_thetaPhiPi0FromOmegaDataHeli, _thetaPhiPi0FromOmegaMcHeli,"_thetaPhiPi0FromOmegaAccCorHeli", "_thetaPhiPi0FromOmegaAccCorHeli");  


  double integralData=_cosOmegaHeliDataHist->Integral();
  DebugMsg << "integralData= " << integralData << endmsg;

  double integralFitted=_cosOmegaHeliFittedHist->Integral();  
  DebugMsg << "integralFittedHists= " << integralFitted << endmsg;

  _cosOmegaHeliFittedHist->Scale(integralData/integralFitted);
  _cosPi0FromOmegaFittedHeli->Scale(integralData/integralFitted);
  _treimanYangFittedHist->Scale(integralData/integralFitted);
  _thetaPhiPi0FromOmegaFittedHeli->Scale(integralData/integralFitted);
  _prodDecThetaPhiPi0FromOmegaFittedHeli->Scale(integralData/integralFitted);

  _cosOmegaFittedAccCorHist=doAccCor(_cosOmegaHeliFittedHist, _cosOmegaHeliMcHist, "_cosOmegaFittedAccCorHist", "cos(#Theta) #omega heli fitted acc cor");
  _cosPi0FromOmegaFittedAccCorHeli=doAccCor(_cosPi0FromOmegaFittedHeli, _cosPi0FromOmegaMcHeli, "_cosPi0FromOmegaFittedAccCorHeli", "cos(#Theta) #pi^0 heli (#omega decay) fitted acc cor");
  _treimanYangFittedAccCorHist=doAccCor(_treimanYangFittedHist, _treimanYangMcHist, "_treimanYangFittedAccCorHist", "_treimanYangFittedAccCorHist");
  _thetaPhiPi0FromOmegaFittedAccCorHeli=doAccCor(_thetaPhiPi0FromOmegaFittedHeli, _thetaPhiPi0FromOmegaMcHeli,"_thetaPhiPi0FromOmegaFittedAccCorHeli", "_thetaPhiPi0FromOmegaFittedAccCorHeli"); 
}

OmegaPiHistLS::~OmegaPiHistLS()
{
  _theTFile->Write();
  _theTFile->Close();
}

void OmegaPiHistLS::initRootStuff(const std::string &thePathToRootFile)
{
  std::stringstream lmaxStrStr;
  std::stringstream pbarMomStrStr;

  lmaxStrStr << _lmax;
  pbarMomStrStr << _pbarmom;

  std::string rootFileName= thePathToRootFile; //"./OmegaPi0Fit_lmax"+lmaxStrStr.str()+"_mom"+pbarMomStrStr.str()+".root";

  _theTFile=new TFile(rootFileName.c_str(),"recreate");
  _cosOmegaHeliDataHist= new TH1F("_cosOmegaHeliDataHist","cos(#Theta) #omega heli data",101, -1.0, 1.0);
  _cosOmegaHeliMcHist= new TH1F("_cosOmegaHeliMcHist","cos(#Theta) #omega heli MC",101, -1.0, 1.0);
  _cosOmegaHeliFittedHist= new TH1F("_cosOmegaHeliFittedHist","cos(#Theta) #omega heli fit result",101, -1.0, 1.0);
  _cosPi0FromOmegaDataHeli= new TH1F("_cosPi0FromOmegaDataHeli","cos(#Theta) #pi^{0} heli (#omega decay) data",101, -1.0, 1.0);
  _cosPi0FromOmegaMcHeli= new TH1F("_cosPi0FromOmegaMcHeli","cos(#Theta) #pi^0 heli (#omega decay) MC",101, -1.0, 1.0);
  _cosPi0FromOmegaFittedHeli= new TH1F("_cosPi0FromOmegaFittedHeli","cos(#Theta) #pi^0 heli (#omega decay) fit result",101, -1.0, 1.0);
  _treimanYangDataHist= new TH1F("_treimanYangDataHist","Treiman Yang angle data",101, -TMath::Pi(), TMath::Pi());
  _treimanYangMcHist= new TH1F("_treimanYangMcHist","Treiman Yang angle MC",101, -TMath::Pi(), TMath::Pi());
  _treimanYangFittedHist= new TH1F("_treimanYangFittedHist","Treiman Yang angle fit result",101, -TMath::Pi(), TMath::Pi()); 
  _cosPi0FromOmegaDataHeli1= new TH1F("_cosPi0FromOmegaDataHeli1","cos(#Theta) #pi^{0} heli (#omega decay) data 1",101, -1.0, 1.0);

  _thetaPhiPi0FromOmegaDataHeli = new TH2F("_thetaPhiPi0FromOmegaDataHeli","_thetaPhiPi0FromOmegaDataHeli",11,-1,1,11,-TMath::Pi(),TMath::Pi());
  _thetaPhiPi0FromOmegaMcHeli = new TH2F("_thetaPhiPi0FromOmegaMcHeli","_thetaPhiPi0FromOmegaMcHeli",11,-1,1,11,-TMath::Pi(),TMath::Pi());
  _thetaPhiPi0FromOmegaFittedHeli = new TH2F("_thetaPhiPi0FromOmegaFittedHeli","_thetaPhiPi0FromOmegaFittedHeli",11,-1,1,11,-TMath::Pi(),TMath::Pi());

  _prodDecThetaPhiPi0FromOmegaDataHeli = new TH3F("__prodDecThetaPhiPi0FromOmegaDataHeli",
						  "__prodDecThetaPhiPi0FromOmegaDataHeli", 8,-1,1,8,-1,1,8,-TMath::Pi(),TMath::Pi());
  _prodDecThetaPhiPi0FromOmegaMcHeli = new TH3F("__prodDecThetaPhiPi0FromOmegaMcHeli",
						"__prodDecThetaPhiPi0FromOmegaMcHeli", 8,-1,1,8,-1,1,8,-TMath::Pi(),TMath::Pi());
  _prodDecThetaPhiPi0FromOmegaFittedHeli = new TH3F("__prodDecThetaPhiPi0FromOmegaFittedHeli",
						    "__prodDecThetaPhiPi0FromOmegaFittedHeli", 8,-1,1,8,-1,1,8,-TMath::Pi(),TMath::Pi());
  
}

TH1F* OmegaPiHistLS::doAccCor(TH1F* dataHist, TH1F* mcHist, const std::string& name, const std::string& title)
{
  TH1F* accCorHist = (TH1F*)(dataHist->Clone());
  accCorHist->Divide(mcHist);
  accCorHist->SetName(name.c_str());
  accCorHist->SetTitle(title.c_str());
  return accCorHist;
}

TH2F* OmegaPiHistLS::doAccCor(TH2F* dataHist, TH2F* mcHist, const std::string& name, const std::string& title)
{
  TH2F* accCorHist = (TH2F*)(dataHist->Clone());
  accCorHist->Divide(mcHist);
  accCorHist->SetName(name.c_str());
  accCorHist->SetTitle(title.c_str());
  return accCorHist;
}

void OmegaPiHistLS::plotCosOmegaHeli(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight)
{
  Vector4<float> omegaHeli4V=theEvtData->omegaHeliCm4Vec;
  theHisto->Fill(omegaHeli4V.CosTheta(), weight);  
}

void OmegaPiHistLS::plotCosPi0FromOmegaHeli(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight)
{
  Vector4<float> piFromOmegaHeli4V=theEvtData->pi0HeliOmega4Vec;
  theHisto->Fill(piFromOmegaHeli4V.CosTheta(), weight);  
}

void OmegaPiHistLS::plotCosPi0FromOmegaHeli1(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight)
{
  theHisto->Fill(theEvtData->cosPi0HeliOmega4Vec, weight);  
}

void OmegaPiHistLS::plotTreimanYang(TH1F* theHisto, const OmPiEvtDataLS* theEvtData, double weight)
{
  Vector4<float> piFromOmegaHeli4V=theEvtData->pi0HeliOmega4Vec;
  theHisto->Fill(piFromOmegaHeli4V.Phi(), weight);  
}

void OmegaPiHistLS::plotThetaPhiPi0FromOmegaHeli(TH2F* theHisto, const OmPiEvtDataLS* theEvtData, double weight)
{
  Vector4<float> piFromOmegaHeli4V=theEvtData->pi0HeliOmega4Vec;
  theHisto->Fill(piFromOmegaHeli4V.CosTheta(), piFromOmegaHeli4V.Phi(), weight);  
}

void OmegaPiHistLS::plotProdDecThetaPhiPi0FromOmegaHeli(TH3F* theHisto, const OmPiEvtDataLS* theEvtData, double weight)
{
  Vector4<float> piFromOmegaHeli4V=theEvtData->pi0HeliOmega4Vec;
  Vector4<float> omegaHeli4V=theEvtData->omegaHeliCm4Vec;
  theHisto->Fill(omegaHeli4V.CosTheta(),piFromOmegaHeli4V.CosTheta(), piFromOmegaHeli4V.Phi(), weight);  
}
