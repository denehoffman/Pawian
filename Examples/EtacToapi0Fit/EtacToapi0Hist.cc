#include <getopt.h>
#include <fstream>
#include <string>
#include "Examples/EtacToapi0Fit/EtacToapi0Hist.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0EventList.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
#include "ErrLogger/ErrLineLog.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

EtacToapi0Hist::EtacToapi0Hist(const EtacToapi0EventList* theEvtList, const fitParamVal& theParamVal) :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _invpietaDataHist(0),
  _invpietaMcHist(0),
  _invpietaFittedHist(0),
  _invpipiDataHist(0),
  _invpipiMcHist(0),
  _invpipiFittedHist(0)
{
  if(0==theEvtList){
    ErrMsg(fatal) <<" EtacToapi0EventList* theEvtList is a 0 pointer !!!!" << endmsg;
  }
  initRootStuff();
  EtacToapi0Lh* theEtacLh=new EtacToapi0Lh(theEvtList);
  const std::vector<evt4Vec> data4Vecs=theEvtList->getDataVecs();

  for (size_t it=0; it<data4Vecs.size(); it++)
    { 
      plotDalitz(_dalitzDataHist, data4Vecs[it], 1.);
      plotInvPiEta(_invpietaDataHist, data4Vecs[it], 1.);
      plotInvPiPi(_invpipiDataHist, data4Vecs[it], 1.);
    } 


  const std::vector<evt4Vec> mc4Vecs=theEvtList->getMcVecs();
  for (size_t it=0; it<mc4Vecs.size(); it++)
    { 
      plotDalitz(_dalitzMcHist, mc4Vecs[it], 1.);
      plotInvPiEta(_invpietaMcHist, mc4Vecs[it], 1.);
      plotInvPiPi(_invpipiMcHist, mc4Vecs[it], 1.);
      
      double evtWeight=theEtacLh->calcEvtIntensity(mc4Vecs[it], theParamVal);
      plotDalitz(_dalitzFittedHist, mc4Vecs[it], evtWeight);
      plotInvPiEta(_invpietaFittedHist, mc4Vecs[it], evtWeight);
      plotInvPiPi(_invpipiFittedHist, mc4Vecs[it], evtWeight);
    }

  //normalize fitted histos

  double integralData=_invpipiDataHist->Integral();
  ErrMsg(debugging) << "integralData= " << integralData << endmsg;
  
  double integralFittedMc=_invpipiFittedHist->Integral();  
  ErrMsg(debugging) << "integralFittedMc= " << integralFittedMc << endmsg;

  double integralMc=_invpipiMcHist->Integral();   
  ErrMsg(debugging) << "integralMc= " << integralMc << endmsg;
  
  _invpipiFittedHist->Scale(integralData/integralFittedMc);
  _invpietaFittedHist->Scale(integralData/integralFittedMc);
  _dalitzFittedHist->Scale(integralData/integralFittedMc);
  _dalitzMcHist->Scale(integralData/integralMc);
  _invpietaMcHist->Scale(integralData/integralMc);
  _invpipiMcHist->Scale(integralData/integralMc); 
  delete theEtacLh;
}

EtacToapi0Hist::~EtacToapi0Hist()
{
  _theTFile->Write();
  _theTFile->Close();
}

void EtacToapi0Hist::initRootStuff()
{
  _theTFile=new TFile("./api0Fit.root","recreate");
  _dalitzDataHist= new TH2F("_dalitzDataHist","dalitzplot data",50, 0., 8.,50, 0., 8.);
  _dalitzMcHist= new TH2F("_dalitzMcHist","dalitzplot mc",50, 0., 8.,50, 0., 8.);
  _dalitzFittedHist= new TH2F("_dalitzFittedHist","dalitzplot fit",50, 0., 8.,50, 0., 8.);
  _invpietaDataHist= new TH1F("_invpietaDataHist","M(#pi #eta) data",120, 0.6, 3.0);
  _invpietaMcHist= new TH1F("_invpietaMcHist","M(#pi #eta) mc",120, 0.6, 3.0);
  _invpietaFittedHist= new TH1F("_invpietaFittedHist","M(#pi #eta) fit",120, 0.6, 3.0);
  _invpipiDataHist= new TH1F("_invpipiDataHist","M(#pi #pi) data",120, 0.2, 2.6);
  _invpipiMcHist= new TH1F("_invpipiMcHist","M(#pi #pi) mc",120, 0.2, 2.6);
  _invpipiFittedHist= new TH1F("_invpipiFittedHist","M(#pi #pi) fit",120, 0.2, 2.6);
}

void EtacToapi0Hist::plotDalitz(TH2F* theHisto, const evt4Vec& theEvt4Vecs, double weight)
{

      Vector4<double> pi0eta4V=theEvt4Vecs.pi0_0_4Vec + theEvt4Vecs.eta_4Vec;
      Vector4<double> pi1eta4V=theEvt4Vecs.pi0_1_4Vec + theEvt4Vecs.eta_4Vec;
      theHisto->Fill(pi0eta4V.M()*pi0eta4V.M(),pi1eta4V.M()*pi1eta4V.M(), weight);
      theHisto->Fill(pi1eta4V.M()*pi1eta4V.M(),pi0eta4V.M()*pi0eta4V.M(), weight);

}

void EtacToapi0Hist::plotInvPiEta (TH1F* theHisto, const evt4Vec& theEvt4Vecs, double weight)
{
      Vector4<double> pi0eta4V=theEvt4Vecs.pi0_0_4Vec + theEvt4Vecs.eta_4Vec;
      Vector4<double> pi1eta4V=theEvt4Vecs.pi0_1_4Vec + theEvt4Vecs.eta_4Vec;
      theHisto->Fill(pi0eta4V.M(), weight);
      theHisto->Fill(pi1eta4V.M(), weight);

}

void EtacToapi0Hist::plotInvPiPi (TH1F* theHisto, const evt4Vec& theEvt4Vecs, double weight)
{

      Vector4<double> pi0pi0=theEvt4Vecs.pi0_0_4Vec + theEvt4Vecs.pi0_1_4Vec;
      theHisto->Fill(pi0pi0.M(), weight);

}
