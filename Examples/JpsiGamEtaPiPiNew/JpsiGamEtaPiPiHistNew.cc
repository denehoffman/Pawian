#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
#include <TStyle.h>

#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiHistNew.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
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


JpsiGamEtaPiPiHistNew::JpsiGamEtaPiPiHistNew(boost::shared_ptr<const EvtDataBaseListNew> theEvtList) :
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
	_massRange(make_pair(0,100))
 {
   if(0==theEvtList){
     Alert <<"theEvtList is a 0 pointer !!!!" ;  // << endmsg;
     exit(1);
   }

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

JpsiGamEtaPiPiHistNew::JpsiGamEtaPiPiHistNew(boost::shared_ptr<AbsLhNew> theJpsiGamEtaPiPiLh, fitParamsNew& theFitParams)
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
	_massRange(make_pair(0,100))
{

  if(0==theJpsiGamEtaPiPiLh){
    Alert <<"JpsiGamEtaPiPiLh* theJpsiGamEtaPiPiLh is a 0 pointer !!!!" ;  // << endmsg;
    exit(1);
  }


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

   double integralData=(double) theEvtList->getDataVecs().size();
   Info <<"No of fit data events  " << integralData << endmsg;   

   double integralFitted=  _EtaPiPiMassFittedHist->Integral();
   Info <<"No of MC events " << theEvtList->getMcVecs().size() << endmsg; 
   Info <<"No of fit events " << integralFitted << endmsg;  

   double scalingFactor=integralData/theEvtList->getMcVecs().size();  
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
	gStyle->SetPalette(1);
	std::stringstream ss;
	ss << _massRange.first << "_" << _massRange.second;

	TCanvas* c_output=new TCanvas( "c_output","c_output",1000,600);
	c_output->Divide(2,2);
	c_output->cd(1);
	_EtaPiMassDataHist->SetLineColor(kRed);
	_EtaPiMassDataHist->DrawCopy("e");
	_EtaPiMassFittedHist->DrawCopy("same");
	c_output->cd(2);
	_PipPimMassDataHist->SetLineColor(kRed);
	_PipPimMassDataHist->DrawCopy("e");
	_PipPimMassFittedHist->DrawCopy("same");
	c_output->cd(3);
	_dalitzDataHist->DrawCopy("colz");
	c_output->cd(4);
	_dalitzFittedHist->DrawCopy("colz");
/*	c_output->cd(3);
	_costEta_EtaPipHeliDataHist->SetLineColor(kRed);
	_costEta_EtaPipHeliDataHist->DrawCopy("e");
	_costEta_EtaPipHeliFittedHist->DrawCopy("same");
	c_output->cd(4);
	_phiEta_EtaPipHeliDataHist->SetLineColor(kRed);
	_phiEta_EtaPipHeliDataHist->DrawCopy("e");
	_phiEta_EtaPipHeliFittedHist->DrawCopy("same");
	c_output->cd(5);
	_costPip_PipPimHeliDataHist->SetLineColor(kRed);
	_costPip_PipPimHeliDataHist->DrawCopy("e");
	_costPip_PipPimHeliFittedHist->DrawCopy("same");
	c_output->cd(6);
	_phiPip_PipPimHeliDataHist->SetLineColor(kRed);
	_phiPip_PipPimHeliDataHist->DrawCopy("e");
	_phiPip_PipPimHeliFittedHist->DrawCopy("same");
*/

	TString outputfilename;
	outputfilename = "JpsiGamEtaPiPi_" + ss.str() + ".pdf";
	c_output->Print(outputfilename);
	c_output->Close();
	_theTFile->Write();
	_theTFile->Close();
}

void JpsiGamEtaPiPiHistNew::initRootStuff(){
	std::string rootFileName="JpsiGamEtaPiPi.root";

	_theTFile=new TFile(rootFileName.c_str(),"RECREATE");

	double xmin;
	double xmax;
//	int xbins = 50;
//	int ybins=xbins;
//	double ymin=xmin;
//	double ymax=xmax;
//	xmin=0.8;
//	xmax=10.;
	_dalitzDataHist= new TH2F("_dalitzDataHist","Dpl #pi^{+}#pi^{-}#eta data", 50, 0.65*0.65, 2.0*2.0, 50, 0.65*0.65, 2.0*2.0 );
	_dalitzMcHist= new TH2F("_dalitzMcHist","Dpl #pi^{+}#pi^{-}#eta MC", 50, 0.65*0.65, 2.0*2.0, 50, 0.65*0.65, 2.0*2.0);
	_dalitzFittedHist= new TH2F("_dalitzFittedHist","Dpl #pi^{+}#pi^{-}#eta fit", 50, 0.65*0.65, 2.0*2.0, 50, 0.65*0.65, 2.0*2.0 );
  
	int nbins=100;
	xmin=0.7;
	xmax=3.2;
	_EtaPiPiMassDataHist   = new TH1F("_EtaPiPiMassDataHist"  , "_EtaPiPiMassDataHist"  , nbins, xmin, xmax);
	_EtaPiPiMassMcHist     = new TH1F("_EtaPiPiMassMcHist"    , "_EtaPiPiMassMcHist"    , nbins, xmin, xmax);
	_EtaPiPiMassFittedHist = new TH1F("_EtaPiPiMassFittedHist", "_EtaPiPiMassFittedHist", nbins, xmin, xmax);
  
	xmin=0.5;
	xmax=2.5;
	_EtaPiMassDataHist   = new TH1F("_EtaPiMassDataHist"  , "_EtaPiMassDataHist"  , nbins, xmin, xmax);
	_EtaPiMassMcHist     = new TH1F("_EtaPiMassMcHist"    , "_EtaPiMassMcHist"    , nbins, xmin, xmax);
	_EtaPiMassFittedHist = new TH1F("_EtaPiMassFittedHist", "_EtaPiMassFittedHist", nbins, xmin, xmax);

	xmin=0.2;
	xmax=2.2;
	_PipPimMassDataHist   = new TH1F("_PipPimMassDataHist"  , "_PipPimMassDataHist"  , nbins, xmin, xmax);
	_PipPimMassMcHist     = new TH1F("_PipPimMassMcHist"    , "_PipPimMassMcHist"    , nbins, xmin, xmax);
	_PipPimMassFittedHist = new TH1F("_PipPimMassFittedHist", "_PipPimMassFittedHist", nbins, xmin, xmax);
  
	_costEta_EtaPipHeliDataHist  = new TH1F("_costEta_EtaPipHeliDataHist"  , "cos(#Theta_{#eta}) #eta#piHeli data", 100, -1., 1.);
	_costEta_EtaPipHeliMcHist    = new TH1F("_costEta_EtaPipHeliMcHist"    , "cos(#Theta_{#eta}) #eta#piHeli Mc"  , 100, -1., 1.);
	_costEta_EtaPipHeliFittedHist= new TH1F("_costEta_EtaPipHeliFittedHist", "cos(#Theta_{#eta}) #eta#piHeli Fit" , 100, -1, 1);
	_phiEta_EtaPipHeliDataHist  = new TH1F("_phiEta_EtaPipHeliDataHist"  , "#Phi_{#eta} #eta#piHeli data", 100, -TMath::Pi(), TMath::Pi());
	_phiEta_EtaPipHeliMcHist    = new TH1F("_phiEta_EtaPipHeliMcHist"    , "#Phi_{#eta} #eta#piHeli Mc"  , 100, -TMath::Pi(), TMath::Pi());
	_phiEta_EtaPipHeliFittedHist= new TH1F("_phiEta_EtaPipHeliFittedHist", "#Phi_{#eta} #eta#piHeli Fit" , 100, -TMath::Pi(), TMath::Pi());

	_costPip_PipPimHeliDataHist  = new TH1F("_costPip_PipPimHeliDataHist"  , "cos(#Theta_{#pi^{+}}) #pi^{+}#pi^{-}Heli data", 100, -1., 1.);
	_costPip_PipPimHeliMcHist    = new TH1F("_costPip_PipPimHeliMcHist"    , "cos(#Theta_{#pi^{+}}) #pi^{+}#pi^{-}Heli Mc"  , 100, -1., 1.);
	_costPip_PipPimHeliFittedHist= new TH1F("_costPip_PipPimHeliFittedHist", "cos(#Theta_{#pi^{+}}) #pi^{+}#pi^{-}Heli Fit" , 100, -1, 1);
	_phiPip_PipPimHeliDataHist  = new TH1F("_phiPip_PipPimHeliDataHist"  , "#Phi_{#pi^{+}} #pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());
	_phiPip_PipPimHeliMcHist    = new TH1F("_phiPip_PipPimHeliMcHist"    , "#Phi_{#pi^{+}} #pi^{+}#pi^{-}Heli Mc"  , 100, -TMath::Pi(), TMath::Pi());
	_phiPip_PipPimHeliFittedHist= new TH1F("_phiPip_PipPimHeliFittedHist", "#Phi_{#pi^{+}} #pi^{+}#pi^{-}Heli Fit" , 100, -TMath::Pi(), TMath::Pi());

	_costPiPi_EtaPiPiHeliDataHist   = new TH1F("_costPiPi_EtaPiPiHeliDataHist"   , "cos(#Theta_{#pi^{+}#pi^{-}}) #eta#pi^{+}#pi^{-}Heli data", 100, -1., 1.);
	_costPiPi_EtaPiPiHeliMcHist     = new TH1F("_costPiPi_EtaPiPiHeliMcHist"     , "cos(#Theta_{#pi^{+}#pi^{-}}) #eta#pi^{+}#pi^{-}Heli Mc", 100, -1., 1.);
	_costPiPi_EtaPiPiHeliFittedHist = new TH1F("_costPiPi_EtaPiPiHeliFittedHist" , "cos(#Theta_{#pi^{+}#pi^{-}}) #eta#pi^{+}#pi^{-}Heli Fit", 100, -1., 1.);
	_phiPiPi_EtaPiPiHeliDataHist   = new TH1F("_phiPiPi_EtaPiPiHeliDataHist"   , "#Phi_{#pi^{+}#pi^{-}} #eta#pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());
	_phiPiPi_EtaPiPiHeliMcHist     = new TH1F("_phiPiPi_EtaPiPiHeliMcHist"     , "#Phi_{#pi^{+}#pi^{-}} #eta#pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());
	_phiPiPi_EtaPiPiHeliFittedHist = new TH1F("_phiPiPi_EtaPiPiHeliFittedHist" , "#Phi_{#pi^{+}#pi^{-}} #eta#pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());

	_costEtaPi_EtaPiPiHeliDataHist   = new TH1F("_costEtaPi_EtaPiPiHeliDataHist"   , "cos(#Theta_{#eta#pi^{+-}}) #eta#pi^{+}#pi^{-}Heli data", 100, -1., 1.);
	_costEtaPi_EtaPiPiHeliMcHist     = new TH1F("_costEtaPi_EtaPiPiHeliMcHist"   , "cos(#Theta_{#eta#pi^{+-}}) #eta#pi^{+}#pi^{-}Heli data", 100, -1., 1.);
	_costEtaPi_EtaPiPiHeliFittedHist = new TH1F("_costEtaPi_EtaPiPiHeliFittedHist"   , "cos(#Theta_{#eta#pi^{+-}}) #eta#pi^{+}#pi^{-}Heli data", 100, -1., 1.);
	_phiEtaPi_EtaPiPiHeliDataHist   = new TH1F("_phiEtaPi_EtaPiPiHeliDataHist"   , "#Phi_{#eta#pi^{+-}} #eta#pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());
	_phiEtaPi_EtaPiPiHeliMcHist     = new TH1F("_phiEtaPi_EtaPiPiHeliMcHist"   , "#Phi_{#eta#pi^{+-}} #eta#pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());
	_phiEtaPi_EtaPiPiHeliFittedHist = new TH1F("_phiEtaPi_EtaPiPiHeliFittedHist"   , "#Phi_{#eta#pi^{+-}} #eta#pi^{+}#pi^{-}Heli data", 100, -TMath::Pi(), TMath::Pi());

	_costGamCmDataHist= new TH1F("_costGamCmDataHist", "cos(#Theta_{#gamma}) CM data", 100, -1., 1.);
	_costGamCmMcHist= new TH1F("_costGamCmMcHist", "cos(#Theta_{#gamma}) CM Mc", 100, -1., 1.);
	_costGamCmFittedHist= new TH1F("_costGamCmFittedHist", "cos(#Theta_{#gamma}) CM Fit", 100, -1, 1);


	std::string tupleVariables = "mEtaPipPim:mEtaPi:mPipPim:EtaPipPimCostTheta:gamCosTheta:EtaPiCosTheta:PipPimCosTheta:PipCosTheta:PimCosTheta:decPlaneChi:testHeli:weight";
  
	_dataTuple=new TNtuple("_dataTuple", "data ntuple", tupleVariables.data());
	_mcTuple=new TNtuple("_mcTuple", "mc ntuple", tupleVariables.data());

	//	_massIndepTuple = new TNtuple("_massIndepTuple","results from mass indep. fits", ("mass:eta:etaErr:f0:f0Err:f1:f1Err:f2:f2Err:eta2:eta2Err:etaReal:etaImg")   );
}

void JpsiGamEtaPiPiHistNew::plotDalitz(TH2F* theHisto,  EvtDataNew* theData, double weight){
	Vector4<double> v4EtaPip = theData->FourVecsDec["EtaPip_HeliPsi"];
	Vector4<double> v4EtaPim = theData->FourVecsDec["EtaPim_HeliPsi"];
	theHisto->Fill( v4EtaPip.M2(), v4EtaPim.M2(), weight);
	theHisto->Fill( v4EtaPim.M2(), v4EtaPip.M2(), weight);
}

void JpsiGamEtaPiPiHistNew::plotEtaPipPim(TH1F* theHisto, EvtDataNew* theData, double weight){
	Vector4<double>& v4 = theData->FourVecsDec["EtaPipPim_HeliPsi"];
	theHisto->Fill( v4.M(), weight );
}

void JpsiGamEtaPiPiHistNew::plotEtaPi(TH1F* theHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["EtaPip_HeliPsi"];
	theHisto->Fill( v4.M(), weight );
	v4 = theData->FourVecsDec["EtaPim_HeliPsi"];
	theHisto->Fill( v4.M(), weight );
}
void JpsiGamEtaPiPiHistNew::plotPipPim(TH1F* theHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["PipPim_HeliPsi"];
	theHisto->Fill( v4.M(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostPhiEta(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["Eta_HeliEtaPip"];
	theCostHisto->Fill( v4.CosTheta(), weight );
	thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostPhiPip(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["Pip_HeliPipPim"];
	theCostHisto->Fill( v4.CosTheta(), weight );
	thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostGam(TH1F* theCostHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["gamma_HeliPsi"];
	theCostHisto->Fill( v4.CosTheta(), weight );
}

void JpsiGamEtaPiPiHistNew::plotCostPhi_PiPiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["PipPim_HeliEtaPipPim"];
	theCostHisto->Fill( v4.CosTheta(), weight);
	thePhiHisto->Fill( v4.Phi(), weight);
}

void JpsiGamEtaPiPiHistNew::plotCostPhi_EtaPiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, EvtDataNew* theData, double weight){
	Vector4<double> v4 = theData->FourVecsDec["EtaPip_HeliEtaPipPim"];
	theCostHisto->Fill( v4.CosTheta(), weight);
	thePhiHisto->Fill( v4.Phi(), weight);
	v4 = theData->FourVecsDec["EtaPim_HeliEtaPipPim"];
	theCostHisto->Fill( v4.CosTheta(), weight);
	thePhiHisto->Fill( v4.Phi(), weight);
}



void  JpsiGamEtaPiPiHistNew::fillTuple( TNtuple* theTuple, EvtDataNew* theData, double weight){
	Vector4<double> V4_EtaPipPim_HeliPsi = theData->FourVecsDec["EtaPipPim_HeliPsi"]  ;
	Vector4<double> V4_EtaPip_HeliPsi= theData->FourVecsDec["EtaPip_HeliPsi"] ;
	Vector4<double> V4_PipPim_HeliPsi= theData->FourVecsDec["PipPim_HeliPsi"] ;
	Vector4<double> V4_gamma_HeliPsi= theData->FourVecsDec["gamma_HeliPsi"] ;
  
	Vector4<double> V4_EtaPip_HeliEtaPipPim= theData->FourVecsDec["EtaPip_HeliEtaPipPim"] ;
	Vector4<double> V4_PipPim_HeliEtaPipPim= theData->FourVecsDec["PipPim_HeliEtaPipPim"] ;
	Vector4<double> V4_Eta_HeliEtaPip= theData->FourVecsDec["Eta_HeliEtaPip"] ;
	Vector4<double> V4_Pip_HeliPipPim= theData->FourVecsDec["Pip_HeliPipPim"] ;
  
//	Vector4<double>& V4_Eta_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliPsi"] ;
//	Vector4<double>& V4_Pip_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPsi"] ;
//	Vector4<double>& V4_Pim_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliPsi"] ;
  
  
	Vector4<double> V4_normPipPimDecHeliEtaPipPim = theData->FourVecsDec["normPipPimDecHeliEtaPipPim"];
	Vector4<double> V4_normEtaPipDecHeliEtaPipPim = theData->FourVecsDec["normEtaPipDecHeliEtaPipPim"];

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





