#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <TStyle.h>

#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiHistNew.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiEventListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"

#include "PwaUtils/KinUtils.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "ErrLogger/ErrLogger.hh"


JpsiGamEtaPiPiHistNew::JpsiGamEtaPiPiHistNew(boost::shared_ptr<const EvtDataBaseListNew> theEvtList, std::pair<double, double> theMassRange) :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _EtaPiPiMassDataHist(0),
  _EtaPiPiMassMcHist(0),
  _EtaPiPiMassFittedHist(0),
  _PipPimMassDataHist(0),
  _PipPimMassMcHist(0),
  _PipPimMassFittedHist(0),
  _EtaPiMassDataHist(0),
  _EtaPiMassMcHist(0),
  _EtaPiMassFittedHist(0),
  _costEta_EtaPipHeliDataHist(0),
  _costEta_EtaPipHeliMcHist(0),
  _costEta_EtaPipHeliFittedHist(0),
  _phiEta_EtaPipHeliDataHist(0),
  _phiEta_EtaPipHeliMcHist(0),
  _phiEta_EtaPipHeliFittedHist(0),
  _costPip_PipPimHeliDataHist(0),
  _costPip_PipPimHeliMcHist(0),
  _costPip_PipPimHeliFittedHist(0),
  _phiPip_PipPimHeliDataHist(0),
  _phiPip_PipPimHeliMcHist(0),
  _phiPip_PipPimHeliFittedHist(0),
  _costPiPi_EtaPiPiHeliDataHist(0),
  _costPiPi_EtaPiPiHeliMcHist(0),
  _costPiPi_EtaPiPiHeliFittedHist(0),
  _phiPiPi_EtaPiPiHeliDataHist(0),
  _phiPiPi_EtaPiPiHeliMcHist(0),
  _phiPiPi_EtaPiPiHeliFittedHist(0),
  _costEtaPi_EtaPiPiHeliDataHist(0),
  _costEtaPi_EtaPiPiHeliMcHist(0),
  _costEtaPi_EtaPiPiHeliFittedHist(0),
  _phiEtaPi_EtaPiPiHeliDataHist(0),
  _phiEtaPi_EtaPiPiHeliMcHist(0),
  _phiEtaPi_EtaPiPiHeliFittedHist(0),
  _costGamCmDataHist(0),
  _costGamCmMcHist(0),
  _costGamCmFittedHist(0),
  _dataTuple(0),
  _mcTuple(0),
  _massIndepTuple(0),
  _massRange(make_pair(0,100)),
  _integralFitted(0),
  _integralData(0),
  _integralMc(0)
{
  if(0==theEvtList){
    Alert <<"theEvtList is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }
  
  _massRange = theMassRange;
  
  initRootStuff();
  
  const std::vector<EvtDataNew*> dataList=theEvtList->getDataVecs();
  
  std::vector<EvtDataNew*>::const_iterator it=dataList.begin();
  
  
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotEtaPipPim(_EtaPiPiMassDataHist, (*it), 1.  );
      plotEtaPi(_EtaPiMassDataHist, (*it), 1. );
      plotPipPim(_PipPimMassDataHist, (*it), 1. );
      plotCostPhiEta(_costEta_EtaPipHeliDataHist , _phiEta_EtaPipHeliDataHist ,(*it), 1. );
      plotCostPhiPip(_costPip_PipPimHeliDataHist, _phiPip_PipPimHeliDataHist,(*it), 1. );
      plotCostGam(_costGamCmDataHist,(*it), 1. );
      plotCostPhiPip(_costPiPi_EtaPiPiHeliDataHist, _phiPiPi_EtaPiPiHeliDataHist,(*it), 1. );
      plotCostPhiEta(_costEtaPi_EtaPiPiHeliDataHist , _phiEtaPi_EtaPiPiHeliDataHist ,(*it), 1. );
      //       fillTuple(_dataTuple, (*it), 1.);
      ++it;
    }
  
  
  const std::vector<EvtDataNew*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    {
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotEtaPipPim(_EtaPiPiMassMcHist, (*it), 1.  );
      plotEtaPi(_EtaPiMassMcHist, (*it), 1. );
      plotPipPim(_PipPimMassMcHist, (*it), 1. );
      plotCostPhiEta(_costEta_EtaPipHeliMcHist , _phiEta_EtaPipHeliMcHist ,(*it), 1. );
      plotCostPhiPip(_costPip_PipPimHeliMcHist, _phiPip_PipPimHeliMcHist,(*it), 1. );
      plotCostGam(_costGamCmMcHist,(*it), 1. );
      plotCostPhiPip(_costPiPi_EtaPiPiHeliMcHist, _phiPiPi_EtaPiPiHeliMcHist,(*it), 1. );
      plotCostPhiEta(_costEtaPi_EtaPiPiHeliMcHist , _phiEtaPi_EtaPiPiHeliMcHist ,(*it), 1. );
      //       fillTuple(_mcTuple, (*it), 1.);
      ++it;
    }
  
  
}

JpsiGamEtaPiPiHistNew::JpsiGamEtaPiPiHistNew(boost::shared_ptr<AbsLhNew> theJpsiGamEtaPiPiLh, fitParamsNew& theFitParams,std::pair<double, double> theMassRange)
  :
  _theTFile(0),
  _dalitzDataHist(0),
  _dalitzMcHist(0),
  _dalitzFittedHist(0),
  _EtaPiPiMassDataHist(0),
  _EtaPiPiMassMcHist(0),
  _EtaPiPiMassFittedHist(0),
  _PipPimMassDataHist(0),
  _PipPimMassMcHist(0),
  _PipPimMassFittedHist(0),
  _EtaPiMassDataHist(0),
  _EtaPiMassMcHist(0),
  _EtaPiMassFittedHist(0),
  _costEta_EtaPipHeliDataHist(0),
  _costEta_EtaPipHeliMcHist(0),
  _costEta_EtaPipHeliFittedHist(0),
  _phiEta_EtaPipHeliDataHist(0),
  _phiEta_EtaPipHeliMcHist(0),
  _phiEta_EtaPipHeliFittedHist(0),
  _costPip_PipPimHeliDataHist(0),
  _costPip_PipPimHeliMcHist(0),
  _costPip_PipPimHeliFittedHist(0),
  _phiPip_PipPimHeliDataHist(0),
  _phiPip_PipPimHeliMcHist(0),
  _phiPip_PipPimHeliFittedHist(0),
  _costPiPi_EtaPiPiHeliDataHist(0),
  _costPiPi_EtaPiPiHeliMcHist(0),
  _costPiPi_EtaPiPiHeliFittedHist(0),
  _phiPiPi_EtaPiPiHeliDataHist(0),
  _phiPiPi_EtaPiPiHeliMcHist(0),
  _phiPiPi_EtaPiPiHeliFittedHist(0),
  _costEtaPi_EtaPiPiHeliDataHist(0),
  _costEtaPi_EtaPiPiHeliMcHist(0),
  _costEtaPi_EtaPiPiHeliFittedHist(0),
  _phiEtaPi_EtaPiPiHeliDataHist(0),
  _phiEtaPi_EtaPiPiHeliMcHist(0),
  _phiEtaPi_EtaPiPiHeliFittedHist(0),
  _costGamCmDataHist(0),
  _costGamCmMcHist(0),
  _costGamCmFittedHist(0),
  _dataTuple(0),
  _mcTuple(0),
  _massIndepTuple(0),
  _massRange(make_pair(0,100)),
  _integralFitted(0),
  _integralData(0),
  _integralMc(0)
{
  
  if(0==theJpsiGamEtaPiPiLh){
    Alert <<"JpsiGamEtaPiPiLh* theJpsiGamEtaPiPiLh is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }
  
  _massRange = theMassRange;
  
  initRootStuff();
  
  boost::shared_ptr<const EvtDataBaseListNew> theEvtList=theJpsiGamEtaPiPiLh->getEventList();
  const std::vector<EvtDataNew*> dataList=theEvtList->getDataVecs();  
  
  std::vector<EvtDataNew*>::const_iterator it=dataList.begin();
  
  
  while(it!=dataList.end())
    {
      plotDalitz(_dalitzDataHist, (*it), 1.);
      plotEtaPipPim(_EtaPiPiMassDataHist, (*it), 1.  );
      plotEtaPi(_EtaPiMassDataHist, (*it), 1. );
      plotPipPim(_PipPimMassDataHist, (*it), 1. );
      plotCostPhiEta(_costEta_EtaPipHeliDataHist , _phiEta_EtaPipHeliDataHist ,(*it), 1. );
      plotCostPhiPip(_costPip_PipPimHeliDataHist, _phiPip_PipPimHeliDataHist,(*it), 1. );
      plotCostGam(_costGamCmDataHist,(*it), 1. );
      plotCostPhiPip(_costPiPi_EtaPiPiHeliDataHist, _phiPiPi_EtaPiPiHeliDataHist,(*it), 1. );
      plotCostPhiEta(_costEtaPi_EtaPiPiHeliDataHist , _phiEtaPi_EtaPiPiHeliDataHist ,(*it), 1. );
      //       fillTuple(_dataTuple, (*it), 1.);
      ++it;
    }
  
  
  const std::vector<EvtDataNew*> mcList=theEvtList->getMcVecs();
  it=mcList.begin();
  while(it!=mcList.end())
    {
      plotDalitz(_dalitzMcHist, (*it), 1.);
      plotEtaPipPim(_EtaPiPiMassMcHist, (*it), 1.  );
      plotEtaPi(_EtaPiMassMcHist, (*it), 1. );
      plotPipPim(_PipPimMassMcHist, (*it), 1. );
      plotCostPhiEta(_costEta_EtaPipHeliMcHist , _phiEta_EtaPipHeliMcHist ,(*it), 1. );
      plotCostPhiPip(_costPip_PipPimHeliMcHist, _phiPip_PipPimHeliMcHist,(*it), 1. );
      plotCostGam(_costGamCmMcHist,(*it), 1. );
      plotCostPhiPip(_costPiPi_EtaPiPiHeliMcHist, _phiPiPi_EtaPiPiHeliMcHist,(*it), 1. );
      plotCostPhiEta(_costEtaPi_EtaPiPiHeliMcHist , _phiEtaPi_EtaPiPiHeliMcHist ,(*it), 1. );
      //       fillTuple(_mcTuple, (*it), 1.);
      
      double evtWeight= theJpsiGamEtaPiPiLh->calcEvtIntensity( (*it), theFitParams );
      plotDalitz(_dalitzFittedHist, (*it), evtWeight);
      plotEtaPipPim(_EtaPiPiMassFittedHist, (*it), evtWeight  );
      plotEtaPi(_EtaPiMassFittedHist, (*it), evtWeight );
      plotPipPim(_PipPimMassFittedHist, (*it), evtWeight );
      plotCostPhiEta(_costEta_EtaPipHeliFittedHist , _phiEta_EtaPipHeliFittedHist ,(*it), evtWeight );
      plotCostPhiPip(_costPip_PipPimHeliFittedHist, _phiPip_PipPimHeliFittedHist,(*it), evtWeight );
      plotCostGam(_costGamCmFittedHist,(*it), evtWeight );
      plotCostPhiPip(_costPiPi_EtaPiPiHeliFittedHist, _phiPiPi_EtaPiPiHeliFittedHist,(*it), evtWeight );
      plotCostPhiEta(_costEtaPi_EtaPiPiHeliFittedHist , _phiEtaPi_EtaPiPiHeliFittedHist ,(*it), evtWeight );
      ++it;
    }
  
  _integralData = (double) theEvtList->getDataVecs().size();
  _integralMc = (double) theEvtList->getMcVecs().size();
  _integralFitted =  _EtaPiPiMassFittedHist->Integral();

  Info <<"No of fit data events " << _integralData << endmsg;   
  Info <<"No of MC events " << _integralMc << endmsg; 
  Info <<"No of fit events " << _integralFitted << endmsg;  
  
  double scalingFactor=_integralData/_integralMc;  
  
  Info <<"scaling factor  " << scalingFactor << endmsg;
  _dalitzFittedHist->Scale(scalingFactor);
  _EtaPiPiMassFittedHist->Scale(scalingFactor);
  _EtaPiMassFittedHist->Scale(scalingFactor);
  _PipPimMassFittedHist->Scale(scalingFactor);
  _costEta_EtaPipHeliFittedHist->Scale(scalingFactor);
  _phiEta_EtaPipHeliFittedHist->Scale(scalingFactor);
  _costPip_PipPimHeliFittedHist->Scale(scalingFactor);
  _phiPip_PipPimHeliFittedHist->Scale(scalingFactor);
  _costGamCmFittedHist->Scale(scalingFactor);
  _costPiPi_EtaPiPiHeliFittedHist->Scale(scalingFactor);
  _phiPiPi_EtaPiPiHeliFittedHist->Scale(scalingFactor);
  _costEtaPi_EtaPiPiHeliFittedHist->Scale(scalingFactor);
  _phiEtaPi_EtaPiPiHeliFittedHist->Scale(scalingFactor);
}

JpsiGamEtaPiPiHistNew::~JpsiGamEtaPiPiHistNew(){

}

void JpsiGamEtaPiPiHistNew::PrintToPDF(std::string suffix){
  gStyle->SetPalette(1);
  std::stringstream ss;
  ss << _massRange.first << "_" << _massRange.second;

  TCanvas* c_output=new TCanvas( "c_output","c_output",1600,900);
  c_output->Divide(3,2);
  c_output->cd(1);
  _EtaPiMassDataHist->SetLineColor(kRed);
  _EtaPiMassDataHist->DrawCopy("e");
  _EtaPiMassFittedHist->DrawCopy("same");
  c_output->cd(2);
  _PipPimMassDataHist->SetLineColor(kRed);
  _PipPimMassDataHist->DrawCopy("e");
  _PipPimMassFittedHist->DrawCopy("same");
  c_output->cd(3);
  _EtaPiPiMassDataHist->SetLineColor(kRed);
  _EtaPiPiMassDataHist->DrawCopy("e");
  _EtaPiPiMassFittedHist->DrawCopy("same");
  c_output->cd(4);
  _costGamCmDataHist->SetLineColor(kRed);
  _costGamCmDataHist->DrawCopy("e");
  _costGamCmFittedHist->DrawCopy("same");
  c_output->cd(5);
  _dalitzDataHist->DrawCopy("colz");
  c_output->cd(6);
  _dalitzFittedHist->DrawCopy("colz");
  TString outputfilename;
  outputfilename = "JpsiGamEtaPiPi_" + ss.str() + suffix.c_str() + ".pdf";
  c_output->Print(outputfilename);
  c_output->Close();

  TCanvas* c_outputangle=new TCanvas( "c_outputangle","c_outputangle",1600,900);
  c_outputangle->Divide(4,2);
  c_outputangle->cd(1);
  _costEta_EtaPipHeliDataHist->SetLineColor(kRed);
  _costEta_EtaPipHeliDataHist->DrawCopy("e");
  _costEta_EtaPipHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(5);
  _phiEta_EtaPipHeliDataHist->SetLineColor(kRed);
  _phiEta_EtaPipHeliDataHist->DrawCopy("e");
  _phiEta_EtaPipHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(2);
  _costEtaPi_EtaPiPiHeliDataHist->SetLineColor(kRed);
  _costEtaPi_EtaPiPiHeliDataHist->DrawCopy("e");
  _costEtaPi_EtaPiPiHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(6);
  _phiEtaPi_EtaPiPiHeliDataHist->SetLineColor(kRed);
  _phiEtaPi_EtaPiPiHeliDataHist->DrawCopy("e");
  _phiEtaPi_EtaPiPiHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(3);
  _costPip_PipPimHeliDataHist->SetLineColor(kRed);
  _costPip_PipPimHeliDataHist->DrawCopy("e");
  _costPip_PipPimHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(7);
  _phiPip_PipPimHeliDataHist->SetLineColor(kRed);
  _phiPip_PipPimHeliDataHist->DrawCopy("e");
  _phiPip_PipPimHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(4);
  _costPiPi_EtaPiPiHeliDataHist->SetLineColor(kRed);
  _costPiPi_EtaPiPiHeliDataHist->DrawCopy("e");
  _costPiPi_EtaPiPiHeliFittedHist->DrawCopy("same");
  c_outputangle->cd(8);
  _phiPiPi_EtaPiPiHeliDataHist->SetLineColor(kRed);
  _phiPiPi_EtaPiPiHeliDataHist->DrawCopy("e");
  _phiPiPi_EtaPiPiHeliFittedHist->DrawCopy("same");
  TString outputanglefilename;
  outputanglefilename = "JpsiGamEtaPiPi_angle_" + ss.str() + suffix.c_str() + ".pdf";
  c_outputangle->Print(outputanglefilename);
  c_outputangle->Close();

  _theTFile->Write();
  _theTFile->Close();
}



void JpsiGamEtaPiPiHistNew::initRootStuff(){
	std::string rootFileName="JpsiGamEtaPiPi.root";

	rootlogon();
	
	_theTFile=new TFile(rootFileName.c_str(),"RECREATE");

	double xmin;
	double xmax;
        double ymin;
        double ymax;

	int nbins=100;
	xmin=_massRange.first-0.05;
	xmax=_massRange.second+0.05;
	_EtaPiPiMassDataHist   = new TH1F("_EtaPiPiMassDataHist"  , "_EtaPiPiMassDataHist"  , nbins, xmin, xmax);
	_EtaPiPiMassMcHist     = new TH1F("_EtaPiPiMassMcHist"    , "_EtaPiPiMassMcHist"    , nbins, xmin, xmax);
	_EtaPiPiMassFittedHist = new TH1F("_EtaPiPiMassFittedHist", "_EtaPiPiMassFittedHist", nbins, xmin, xmax);
  
	nbins=50;
	xmin=0.65;
	xmax=_massRange.second-0.1;
	ymin=0.65;
        ymax=_massRange.second-0.1;

	_EtaPiMassDataHist   = new TH1F("_EtaPiMassDataHist"  , "_EtaPiMassDataHist"  , nbins, xmin, xmax);
	_EtaPiMassMcHist     = new TH1F("_EtaPiMassMcHist"    , "_EtaPiMassMcHist"    , nbins, xmin, xmax);
	_EtaPiMassFittedHist = new TH1F("_EtaPiMassFittedHist", "_EtaPiMassFittedHist", nbins, xmin, xmax);

        _dalitzDataHist   = new TH2F("_dalitzDataHist"  ,"Dpl #pi^{+}#pi^{-}#eta data; #m^{2}(#eta#pi);#m^{2}(#eta#pi)", nbins, xmin*xmin, xmax*xmax, nbins, ymin*ymin, ymax*ymax );
        _dalitzMcHist     = new TH2F("_dalitzMcHist"    ,"Dpl #pi^{+}#pi^{-}#eta MC  ; #m^{2}(#eta#pi);#m^{2}(#eta#pi)", nbins, xmin*xmin, xmax*xmax, nbins, ymin*ymin, ymax*ymax );
        _dalitzFittedHist = new TH2F("_dalitzFittedHist","Dpl #pi^{+}#pi^{-}#eta fit ; #m^{2}(#eta#pi);#m^{2}(#eta#pi)", nbins, xmin*xmin, xmax*xmax, nbins, ymin*ymin, ymax*ymax );

	xmin=0.25;
	xmax=_massRange.second-0.5;
	_PipPimMassDataHist   = new TH1F("_PipPimMassDataHist"  , "_PipPimMassDataHist"  , nbins, xmin, xmax);
	_PipPimMassMcHist     = new TH1F("_PipPimMassMcHist"    , "_PipPimMassMcHist"    , nbins, xmin, xmax);
	_PipPimMassFittedHist = new TH1F("_PipPimMassFittedHist", "_PipPimMassFittedHist", nbins, xmin, xmax);
  
	_costEta_EtaPipHeliDataHist  = new TH1F("_costEta_EtaPipHeliDataHist"  , "cos(#Theta_{#eta}) #eta#piHeli data", nbins, -1., 1.);
	_costEta_EtaPipHeliMcHist    = new TH1F("_costEta_EtaPipHeliMcHist"    , "cos(#Theta_{#eta}) #eta#piHeli Mc"  , nbins, -1., 1.);
	_costEta_EtaPipHeliFittedHist= new TH1F("_costEta_EtaPipHeliFittedHist", "cos(#Theta_{#eta}) #eta#piHeli Fit" , nbins, -1., 1.);
	_phiEta_EtaPipHeliDataHist   = new TH1F("_phiEta_EtaPipHeliDataHist"   , "#Phi_{#eta} #eta#piHeli data", nbins, -TMath::Pi(), TMath::Pi());
	_phiEta_EtaPipHeliMcHist     = new TH1F("_phiEta_EtaPipHeliMcHist"     , "#Phi_{#eta} #eta#piHeli Mc"  , nbins, -TMath::Pi(), TMath::Pi());
	_phiEta_EtaPipHeliFittedHist = new TH1F("_phiEta_EtaPipHeliFittedHist" , "#Phi_{#eta} #eta#piHeli Fit" , nbins, -TMath::Pi(), TMath::Pi());

	_costPip_PipPimHeliDataHist  = new TH1F("_costPip_PipPimHeliDataHist"  , "cos(#Theta_{#pi^{+}}) #pi^{+}#pi^{-}Heli data", nbins, -1., 1.);
	_costPip_PipPimHeliMcHist    = new TH1F("_costPip_PipPimHeliMcHist"    , "cos(#Theta_{#pi^{+}}) #pi^{+}#pi^{-}Heli Mc"  , nbins, -1., 1.);
	_costPip_PipPimHeliFittedHist= new TH1F("_costPip_PipPimHeliFittedHist", "cos(#Theta_{#pi^{+}}) #pi^{+}#pi^{-}Heli Fit" , nbins, -1., 1.);
	_phiPip_PipPimHeliDataHist   = new TH1F("_phiPip_PipPimHeliDataHist"   , "#Phi_{#pi^{+}} #pi^{+}#pi^{-}Heli data", nbins, -TMath::Pi(), TMath::Pi());
	_phiPip_PipPimHeliMcHist     = new TH1F("_phiPip_PipPimHeliMcHist"     , "#Phi_{#pi^{+}} #pi^{+}#pi^{-}Heli Mc"  , nbins, -TMath::Pi(), TMath::Pi());
	_phiPip_PipPimHeliFittedHist = new TH1F("_phiPip_PipPimHeliFittedHist" , "#Phi_{#pi^{+}} #pi^{+}#pi^{-}Heli Fit" , nbins, -TMath::Pi(), TMath::Pi());

	_costPiPi_EtaPiPiHeliDataHist   = new TH1F("_costPiPi_EtaPiPiHeliDataHist"   , "cos(#Theta_{#pi^{+}#pi^{-}}) #eta#pi^{+}#pi^{-}Heli data", nbins, -1., 1.);
	_costPiPi_EtaPiPiHeliMcHist     = new TH1F("_costPiPi_EtaPiPiHeliMcHist"     , "cos(#Theta_{#pi^{+}#pi^{-}}) #eta#pi^{+}#pi^{-}Heli Mc"  , nbins, -1., 1.);
	_costPiPi_EtaPiPiHeliFittedHist = new TH1F("_costPiPi_EtaPiPiHeliFittedHist" , "cos(#Theta_{#pi^{+}#pi^{-}}) #eta#pi^{+}#pi^{-}Heli Fit" , nbins, -1., 1.);
	_phiPiPi_EtaPiPiHeliDataHist   = new TH1F("_phiPiPi_EtaPiPiHeliDataHist"   , "#Phi_{#pi^{+}#pi^{-}} #eta#pi^{+}#pi^{-}Heli data", nbins, -TMath::Pi(), TMath::Pi());
	_phiPiPi_EtaPiPiHeliMcHist     = new TH1F("_phiPiPi_EtaPiPiHeliMcHist"     , "#Phi_{#pi^{+}#pi^{-}} #eta#pi^{+}#pi^{-}Heli Mc"  , nbins, -TMath::Pi(), TMath::Pi());
	_phiPiPi_EtaPiPiHeliFittedHist = new TH1F("_phiPiPi_EtaPiPiHeliFittedHist" , "#Phi_{#pi^{+}#pi^{-}} #eta#pi^{+}#pi^{-}Heli Fit" , nbins, -TMath::Pi(), TMath::Pi());

	_costEtaPi_EtaPiPiHeliDataHist   = new TH1F("_costEtaPi_EtaPiPiHeliDataHist"  , "cos(#Theta_{#eta#pi^{+-}}) #eta#pi^{+}#pi^{-}Heli data", nbins, -1., 1.);
	_costEtaPi_EtaPiPiHeliMcHist     = new TH1F("_costEtaPi_EtaPiPiHeliMcHist"    , "cos(#Theta_{#eta#pi^{+-}}) #eta#pi^{+}#pi^{-}Heli Mc  ", nbins, -1., 1.);
	_costEtaPi_EtaPiPiHeliFittedHist = new TH1F("_costEtaPi_EtaPiPiHeliFittedHist", "cos(#Theta_{#eta#pi^{+-}}) #eta#pi^{+}#pi^{-}Heli Fit ", nbins, -1., 1.);
	_phiEtaPi_EtaPiPiHeliDataHist   = new TH1F("_phiEtaPi_EtaPiPiHeliDataHist"    , "#Phi_{#eta#pi^{+-}} #eta#pi^{+}#pi^{-}Heli data", nbins, -TMath::Pi(), TMath::Pi());
	_phiEtaPi_EtaPiPiHeliMcHist     = new TH1F("_phiEtaPi_EtaPiPiHeliMcHist"      , "#Phi_{#eta#pi^{+-}} #eta#pi^{+}#pi^{-}Heli Mc"  , nbins, -TMath::Pi(), TMath::Pi());
	_phiEtaPi_EtaPiPiHeliFittedHist = new TH1F("_phiEtaPi_EtaPiPiHeliFittedHist"  , "#Phi_{#eta#pi^{+-}} #eta#pi^{+}#pi^{-}Heli Fit" , nbins, -TMath::Pi(), TMath::Pi());

	nbins = 100;
	_costGamCmDataHist   = new TH1F("_costGamCmDataHist"  , "cos(#Theta_{#gamma}) CM data", nbins, -1., 1.);
	_costGamCmMcHist     = new TH1F("_costGamCmMcHist"    , "cos(#Theta_{#gamma}) CM Mc"  , nbins, -1., 1.);
	_costGamCmFittedHist = new TH1F("_costGamCmFittedHist", "cos(#Theta_{#gamma}) CM Fit" , nbins, -1., 1.);


	std::string tupleVariables = "mEtaPipPim:mEtaPi:mPipPim:EtaPipPimCostTheta:gamCosTheta:EtaPiCosTheta:PipPimCosTheta:PipCosTheta:PimCosTheta:decPlaneChi:testHeli:weight";
  
	_dataTuple=new TNtuple("_dataTuple", "data ntuple", tupleVariables.data());
	_mcTuple=new TNtuple("_mcTuple", "mc ntuple", tupleVariables.data());

	//	_massIndepTuple = new TNtuple("_massIndepTuple","results from mass indep. fits", ("mass:eta:etaErr:f0:f0Err:f1:f1Err:f2:f2Err:eta2:eta2Err:etaReal:etaImg")   );
}

void JpsiGamEtaPiPiHistNew::rootlogon(){
  // Canvas style

  gStyle->SetPadBorderMode(1);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadColor(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetCanvasBorderMode(1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  
  // Fit and stat. format 

  gStyle->SetOptFit(1011);
  gStyle->SetFitFormat("2.5f");
  gStyle->SetOptStat(1000010) ;
  gStyle->SetStatColor(0);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.97);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  
  // Histogram style 

  gStyle->SetTitleColor(1);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleY(0.975);

  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.8,"y");
  gStyle->SetTitleOffset(1.7,"z");

  gStyle->SetLabelOffset(.02,"x");
  gStyle->SetLabelOffset(.025,"y");
  gStyle->SetLabelOffset(.02,"z");

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelFont(42,"z");

  gStyle->SetTickLength(-0.025,"x");
  gStyle->SetTickLength(-0.025,"y");
  gStyle->SetTickLength(-0.025,"z");

  gStyle->SetLabelSize(.0425,"x");
  gStyle->SetLabelSize(.0425,"y");
  gStyle->SetLabelSize(.0425,"z");

  gStyle->SetTitleSize(.0425,"x");
  gStyle->SetTitleSize(.0425,"y");
  gStyle->SetTitleSize(.0425,"z");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.23);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);

  gStyle->SetTitleFont(42);
  gStyle->SetStatFont(42);

  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleFont(42,"z");

  gStyle->SetPalette(1);

  gStyle->SetHistFillStyle(1001);
  gStyle->SetHistLineWidth(1);

}

void JpsiGamEtaPiPiHistNew::plotDalitz(TH2F* theHisto,  EvtDataNew* theData, double weight){
  Vector4<double> v4EtaPip = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi];
	Vector4<double> v4EtaPim = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi];
	theHisto->Fill( v4EtaPip.M2(), v4EtaPim.M2(), weight);
	theHisto->Fill( v4EtaPim.M2(), v4EtaPip.M2(), weight);
}

void JpsiGamEtaPiPiHistNew::plotEtaPipPim(TH1F* theHisto, EvtDataNew* theData, double weight){
	Vector4<double>& v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPipPim_HeliPsi];
	theHisto->Fill( v4.M(), weight );
}

void JpsiGamEtaPiPiHistNew::plotEtaPi(TH1F* theHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi];
	theHisto->Fill( v4.M(), weight );
	v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi];
	theHisto->Fill( v4.M(), weight );
}
void JpsiGamEtaPiPiHistNew::plotPipPim(TH1F* theHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliPsi];
	theHisto->Fill( v4.M(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostPhiEta(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPip];
	theCostHisto->Fill( v4.CosTheta(), weight );
	thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostPhiPip(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pip_HeliPipPim];
	theCostHisto->Fill( v4.CosTheta(), weight );
	thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostGam(TH1F* theCostHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::gamma_HeliPsi];
	theCostHisto->Fill( v4.CosTheta(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostPhi_PiPiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliEtaPipPim];
	theCostHisto->Fill( v4.CosTheta(), weight);
	thePhiHisto->Fill( v4.Phi(), weight);
}

void JpsiGamEtaPiPiHistNew::plotCostPhi_EtaPiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliEtaPipPim];
	theCostHisto->Fill( v4.CosTheta(), weight);
	thePhiHisto->Fill( v4.Phi(), weight);
	v4 = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliEtaPipPim];
	theCostHisto->Fill( v4.CosTheta(), weight);
	thePhiHisto->Fill( v4.Phi(), weight);
}



void  JpsiGamEtaPiPiHistNew::fillTuple( TNtuple* theTuple, EvtDataNew* theData, double weight){
	Vector4<double> V4_EtaPipPim_HeliPsi = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPipPim_HeliPsi]  ;
	Vector4<double> V4_EtaPip_HeliPsi= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi] ;
	Vector4<double> V4_PipPim_HeliPsi= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliPsi] ;
	Vector4<double> V4_gamma_HeliPsi= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::gamma_HeliPsi] ;
  
	Vector4<double> V4_EtaPip_HeliEtaPipPim= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliEtaPipPim] ;
	Vector4<double> V4_PipPim_HeliEtaPipPim= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliEtaPipPim] ;
	Vector4<double> V4_Eta_HeliEtaPip= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPip] ;
	Vector4<double> V4_Pip_HeliPipPim= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::Pip_HeliPipPim] ;
  
//	Vector4<double>& V4_Eta_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliPsi"] ;
//	Vector4<double>& V4_Pip_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPsi"] ;
//	Vector4<double>& V4_Pim_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliPsi"] ;
  
  
	Vector4<double> V4_normPipPimDecHeliEtaPipPim = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::normPipPimDecHeliEtaPipPim];
	Vector4<double> V4_normEtaPipDecHeliEtaPipPim = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::normEtaPipDecHeliEtaPipPim];

	double cosChi=(V4_normPipPimDecHeliEtaPipPim.Px()*V4_normEtaPipDecHeliEtaPipPim.Px()
				+V4_normPipPimDecHeliEtaPipPim.Py()*V4_normEtaPipDecHeliEtaPipPim.Py()
				+V4_normPipPimDecHeliEtaPipPim.Pz()*V4_normEtaPipDecHeliEtaPipPim.Pz())
				/ (V4_normPipPimDecHeliEtaPipPim.P()*V4_normEtaPipDecHeliEtaPipPim.P());

	double chi=acos(fabs(cosChi));
	double thePhiPhiDecayPlaneAngle = chi/TMath::Pi()*180.;
  
	//double thePhiPhiDecayPlaneAngle = decayAngleChi( V4_KsKlKpKm_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi, V4_Ks_HeliPsi, V4_Kl_HeliPsi   );
	//double testHeli = costDecHeli( V4_KsKlKpKm_HeliPsi+V4_gamma_HeliPsi, V4_Ks_HeliPsi+V4_Kl_HeliPsi+V4_Km_HeliPsi+V4_Kp_HeliPsi, V4_Ks_HeliPsi+V4_Kl_HeliPsi );
  
	cout << "FillTuple" << endl;
	theTuple->Fill(
			V4_EtaPipPim_HeliPsi.M(),
			V4_EtaPip_HeliPsi.M(),
			V4_PipPim_HeliPsi.M(),
			V4_EtaPipPim_HeliPsi.CosTheta(),
			V4_gamma_HeliPsi.CosTheta(),
			V4_EtaPip_HeliEtaPipPim.CosTheta(),
			V4_PipPim_HeliEtaPipPim.CosTheta(),
			V4_Eta_HeliEtaPip.CosTheta(),
			V4_Pip_HeliPipPim.CosTheta(),
			thePhiPhiDecayPlaneAngle,
//			testHeli,
			weight
			);
}





