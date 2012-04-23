#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiHist.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiEventList.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiProdLh.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiFitParams.hh"

#include "PwaUtils/KinUtils.hh"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "ErrLogger/ErrLogger.hh"


JpsiGamEtaPiPiHist::JpsiGamEtaPiPiHist(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
	_theTFile(0),
	//_dalitzDataHist(0),
	//_dalitzMcHist(0),
	//_dalitzFittedHist(0),
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
	_costGamCmDataHist(0),
	_costGamCmMcHist(0),
	_costGamCmFittedHist(0),
	_chiDataHist(0),
	_chiMcHist(0),
	_chiFittedHist(0),
	_dataTuple(0),
	_mcTuple(0),
	_massIndepTuple(0),
	_massRange(make_pair(0,100))
 {
   if(0==theEvtList){
     Alert <<"JpsiGamEtaPiPiEventList* theEvtList is a 0 pointer !!!!" ;  // << endmsg;
     exit(1);
   }

   initRootStuff();
  
   const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

   std::vector<EvtData*>::const_iterator it=dataList.begin();


   while(it!=dataList.end())
     {
	   //       plotDalitz(_dalitzDataHist, (*it), 1.);

       plotEtaPipPim(_EtaPiPiMassDataHist, (*it), 1.  );
       plotEtaPi(_EtaPiMassDataHist, (*it), 1. );
       plotPipPim(_PipPimMassDataHist, (*it), 1. );
       plotCostPhiEta(_costEta_EtaPipHeliDataHist , _phiEta_EtaPipHeliDataHist ,(*it), 1. );
       plotCostPhiPip(_costPip_PipPimHeliDataHist, _phiPip_PipPimHeliDataHist,(*it), 1. );
       plotCostGam(_costGamCmDataHist,(*it), 1. );
//       plotCostPhi_PhiPhiHeli(_costPhi_KpKmDataHist, _phiPhi_KpKmDataHist, (*it)->FourVecs[enumJpsiGamEtaPiPiData::V4_KpKm_HeliKsKlKpKm], 1.);
//       fillTuple(_dataTuple, (*it), 1.);
       ++it;
     }

  
   const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
   it=mcList.begin();
   while(it!=mcList.end())
     {
	   //       plotDalitz(_dalitzMcHist, (*it), 1.);
       plotEtaPipPim(_EtaPiPiMassMcHist, (*it), 1.  );
       plotEtaPi(_EtaPiMassMcHist, (*it), 1. );
       plotPipPim(_PipPimMassMcHist, (*it), 1. );
       plotCostPhiEta(_costEta_EtaPipHeliMcHist , _phiEta_EtaPipHeliMcHist ,(*it), 1. );
       plotCostPhiPip(_costPip_PipPimHeliMcHist, _phiPip_PipPimHeliMcHist,(*it), 1. );
       plotCostGam(_costGamCmMcHist,(*it), 1. );
//       plotCostPhi_PhiPhiHeli(_costPhi_KpKmMcHist, _phiPhi_KpKmMcHist, (*it)->FourVecs[enumJpsiGamEtaPiPiMc::V4_KpKm_HeliKsKlKpKm], 1.);
//       fillTuple(_mcTuple, (*it), 1.);
       ++it;
     }
   

 }


JpsiGamEtaPiPiHist::JpsiGamEtaPiPiHist(JpsiGamEtaPiPiProdLh* theJpsiGamEtaPiPiLh, fitParams& fitParam, FitParamErrorMatrix* theErrorMatrix) :
		_theTFile(0),
		//_dalitzDataHist(0),
		//_dalitzMcHist(0),
		//_dalitzFittedHist(0),
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
		_costGamCmDataHist(0),
		_costGamCmMcHist(0),
		_costGamCmFittedHist(0),
		_chiDataHist(0),
		_chiMcHist(0),
		_chiFittedHist(0),
		_dataTuple(0),
		_mcTuple(0),
		_massIndepTuple(0),
		_massRange(make_pair(0,100)){

	if(0==theJpsiGamEtaPiPiLh){
		Alert <<"JpsiGamEtaPiPiLh* theJpsiGamEtaPiPiLh is a 0 pointer !!!!" ;  // << endmsg;
		exit(1);
	}

	initRootStuff();
	_theJpsiGamEtaPiPiLh = theJpsiGamEtaPiPiLh;
	_fitParam = fitParam;
	std::vector<double> data;
	_errMatrix = theErrorMatrix;
}
  

void JpsiGamEtaPiPiHist::fill(){
	boost::shared_ptr<const EvtDataBaseList> theEvtList=_theJpsiGamEtaPiPiLh->getEventList();
	const std::vector<EvtData*> dataList=theEvtList->getDataVecs();

	std::vector<EvtData*>::const_iterator it=dataList.begin();
	while(it!=dataList.end()){
		//		plotDalitz(_dalitzDataHist, (*it), 1.);
		plotEtaPipPim(_EtaPiPiMassDataHist, (*it), 1.  );
		plotEtaPi(_EtaPiMassDataHist, (*it), 1. );
		plotPipPim(_PipPimMassDataHist, (*it), 1. );
		plotCostPhiEta(_costEta_EtaPipHeliDataHist , _phiEta_EtaPipHeliDataHist ,(*it), 1. );
		plotCostPhiPip(_costPip_PipPimHeliDataHist, _phiPip_PipPimHeliDataHist,(*it), 1. );
		plotCostGam(_costGamCmDataHist,(*it), 1. );
		//plotCostPhi_PhiPhiHeli(_costPhi_KpKmDataHist, _phiPhi_KpKmDataHist, (*it)->FourVecs[enumJpsiGamEtaPiPiData::V4_KpKm_HeliKsKlKpKm], 1.);
//		fillTuple(_dataTuple, (*it), 1.);

		++it;
	}

	const std::vector<EvtData*> mcList=theEvtList->getMcVecs();
	it=mcList.begin();
	while(it!=mcList.end()){
		//		plotDalitz(_dalitzMcHist, (*it), 1.);
		plotEtaPipPim(_EtaPiPiMassMcHist, (*it), 1.  );
		plotEtaPi(_EtaPiMassMcHist, (*it), 1. );
		plotPipPim(_PipPimMassMcHist, (*it), 1. );
		plotCostPhiEta(_costEta_EtaPipHeliMcHist , _phiEta_EtaPipHeliMcHist ,(*it), 1. );
		plotCostPhiPip(_costPip_PipPimHeliMcHist, _phiPip_PipPimHeliMcHist,(*it), 1. );
		plotCostGam(_costGamCmMcHist,(*it), 1. );
		//plotCostPhi_PhiPhiHeli(_costPhi_KpKmMcHist, _phiPhi_KpKmMcHist, (*it)->FourVecs[enumJpsiGamEtaPiPiMc::V4_KpKm_HeliKsKlKpKm], 1.);
//		fillTuple(_mcTuple, (*it), 1.);

		double evtWeight= _theJpsiGamEtaPiPiLh->calcEvtIntensity((*it), _fitParam);
		//		plotDalitz(_dalitzFittedHist, (*it), evtWeight);
		plotEtaPipPim(_EtaPiPiMassFittedHist, (*it), evtWeight  );
		plotEtaPi(_EtaPiMassFittedHist, (*it), evtWeight );
		plotPipPim(_PipPimMassFittedHist, (*it), evtWeight );
		plotCostPhiEta(_costEta_EtaPipHeliFittedHist , _phiEta_EtaPipHeliFittedHist ,(*it), evtWeight );
		plotCostPhiPip(_costPip_PipPimHeliFittedHist, _phiPip_PipPimHeliFittedHist,(*it), evtWeight );
		plotCostGam(_costGamCmFittedHist,(*it), evtWeight );
		//plotCostPhi_PhiPhiHeli(_costPhi_KpKmFittedHist, _phiPhi_KpKmFittedHist, (*it)->FourVecs[enumJpsiGamEtaPiPiFitted::V4_KpKm_HeliKsKlKpKm], evtWeight);
		plotChi(_chiFittedHist, (*it), evtWeight  );

//		fillTuple(_mcTuple, (*it), evtWeight);

		++it;
    }


	double integralData=(double) theEvtList->getDataVecs().size();
	Info <<"No of fit data events  " << integralData ;  // << endmsg;
  
	double integralFitted=(double) theEvtList->getMcVecs().size();
	Info <<"No of fit events " << integralFitted ;  // << endmsg;

	Info <<"scaling factor  " << integralData/integralFitted ;  // << endmsg;

	double scaleFactor = integralData/integralFitted;
  
	//_dalitzFittedHist->Scale(scaleFactor);
	_EtaPiPiMassFittedHist->Scale(scaleFactor);
	_EtaPiMassFittedHist->Scale(scaleFactor);
	_PipPimMassFittedHist->Scale(scaleFactor);
	_costEta_EtaPipHeliFittedHist->Scale(scaleFactor);
	_phiEta_EtaPipHeliFittedHist->Scale(scaleFactor);
	_costPip_PipPimHeliFittedHist->Scale(scaleFactor);
	_phiPip_PipPimHeliFittedHist->Scale(scaleFactor);
	_costGamCmFittedHist->Scale(scaleFactor);
	//_costPhi_KpKmFittedHist->Scale(scaleFactor);
	//_phiPhi_KpKmFittedHist->Scale(scaleFactor);
	_chiFittedHist->Scale(scaleFactor);

	/*
	double iEta(0.), iEtaErr(0.), iF0(0.), iF0Err(0.), iEta2(0.), iEta2Err(0.), iF1(0.), iF1Err(0.), iF2(0.), iF2Err(0.);
	double etaReal(0.), etaImg(0.);


	it=mcList.begin();
	while(it!=mcList.end()){
		std::pair<double, double> intensityEvent = make_pair(0.,0.);
		_theJpsiGamEtaPiPiLh->calcComponentIntensity((*it), _fitParam, *_errMatrix, "etacHyp", intensityEvent);
		iEta+= intensityEvent.first*scaleFactor;
		iEtaErr+= intensityEvent.second*scaleFactor;

		_theJpsiGamEtaPiPiLh->calcComponentIntensity((*it), _fitParam, *_errMatrix, "f02020Hyp",intensityEvent);
		iF0 += intensityEvent.first*scaleFactor;
		iF0Err += intensityEvent.second*scaleFactor;

		_theJpsiGamEtaPiPiLh->calcComponentIntensity((*it), _fitParam, *_errMatrix, "f1Hyp", intensityEvent);
		iF1 += intensityEvent.first*scaleFactor;
		iF1Err += intensityEvent.second*scaleFactor;

		_theJpsiGamEtaPiPiLh->calcComponentIntensity((*it), _fitParam, *_errMatrix, "f22300Hyp",intensityEvent);
		iF2 += intensityEvent.first*scaleFactor;
		iF2Err += intensityEvent.second*scaleFactor;

		_theJpsiGamEtaPiPiLh->calcComponentIntensity((*it), _fitParam, *_errMatrix, "eta21870Hyp",intensityEvent);
		iEta2+= intensityEvent.first*scaleFactor;
		iEta2Err += intensityEvent.second*scaleFactor;

		it++;
	}
	double meanMassRange = _massRange.first + 0.5*(_massRange.second-_massRange.first);

	Info << "Events for component eta : " << iEta << " +/- " << iEtaErr ;
	Info << "Events for component f0:   " << iF0 << " +/- " << iF0Err ;
	Info << "Events for component f1:   " << iF1 << " +/- " << iF1Err ;
	Info << "Events for component f2:   " << iF2 << " +/- " << iF2Err ;
	Info << "Events for component eta2: " << iEta2 << " +/- " << iEta2Err ;

	std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamPhi=_fitParam.Phis[paramEnumJpsiGamEtaPiPi::PsiToEtacGamma];
	std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamMag=_fitParam.Mags[paramEnumJpsiGamEtaPiPi::PsiToEtacGamma];
	std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator phiIter = PsiToEtacGamPhi.begin();
	for(; phiIter!=PsiToEtacGamPhi.end(); phiIter++){
		double thePhase = phiIter->second;
		double theMag = PsiToEtacGamMag[phiIter->first];

		etaImg = theMag*cos(thePhase) ;
		etaReal = theMag*sin(thePhase);
	}

	_massIndepTuple->Fill(meanMassRange, iEta, iEtaErr, iF0, iF0Err, iF1, iF1Err, iF2, iF2Err, iEta2, iEta2Err, etaImg, etaReal );
*/
}

JpsiGamEtaPiPiHist::~JpsiGamEtaPiPiHist(){
	_theTFile->Write();
	_theTFile->Close();
}

void JpsiGamEtaPiPiHist::initRootStuff(){
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
//	_dalitzDataHist= new TH2F("_dalitzDataHist","Dpl K+K- K+#pi^{0} data",xbins, xmin, xmax, ybins, ymin, ymax );
//	_dalitzMcHist= new TH2F("_dalitzMcHist","Dpl K+K- K+#pi^{0} MC",xbins, xmin, xmax, ybins, ymin, ymax);
//	_dalitzFittedHist= new TH2F("_dalitzFittedHist","Dpl K+K- K+#pi^{0} fit",xbins, xmin, xmax, ybins, ymin, ymax );

  
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

	_costGamCmDataHist= new TH1F("_costGamCmDataHist", "cos(#Theta_{#gamma}) CM data", 100, -1., 1.);
	_costGamCmMcHist= new TH1F("_costGamCmMcHist", "cos(#Theta_{#gamma}) CM Mc", 100, -1., 1.);
	_costGamCmFittedHist= new TH1F("_costGamCmFittedHist", "cos(#Theta_{#gamma}) CM Fit", 100, -1, 1);

	//	_costPhi_KpKmDataHist= new TH1F("_costPhi_KpKmDataHist", "cos(#Theta_{#phi}) K+ K- data", 100, -1., 1.);
	//	_costPhi_KpKmMcHist= new TH1F("_costPhi_KpKmMcHist", "cos(#Theta_{#phi}) K+ K- Mc", 100, -1., 1.);
	//	_costPhi_KpKmFittedHist= new TH1F("_costPhi_KpKmFittedHist", "cos(#Theta_{#phi}) K+ K- Fit", 100, -1., 1.);

	//	_phiPhi_KpKmDataHist= new TH1F("_phiPhi_KpKmDataHist", "cos(#Phi_{#phi}) K+ K- data", 100, -TMath::Pi(), TMath::Pi());
	//	_phiPhi_KpKmMcHist= new TH1F("_phiPhi_KpKmMcHist", "cos(#Phi_{#phi}) K+ K- Mc", 100, -TMath::Pi(), TMath::Pi());
	//	_phiPhi_KpKmFittedHist= new TH1F("_phiPhi_KpKmFittedHist", "cos(#Phi_{#phi}) K+ K- Fit", 100, -TMath::Pi(), TMath::Pi());

	_chiDataHist= new TH1F("_chiDataHist", "#chi data", 90, 0., 90.);
	_chiMcHist= new TH1F("_chiMcHist", "#chi Mc", 90, 0., 90.);
	_chiFittedHist= new TH1F("_chiFittedHist", "#chi Fit", 90, 0., 90.);

	std::string tupleVariables = "mEtaPipPim:mEtaPi:mPipPim:EtaPipPimCostTheta:gamCosTheta:EtaPiCosTheta:PipPimCosTheta:PipCosTheta:PimCosTheta:decPlaneChi:testHeli:weight";
  
	_dataTuple=new TNtuple("_dataTuple", "data ntuple", tupleVariables.data());
	_mcTuple=new TNtuple("_mcTuple", "mc ntuple", tupleVariables.data());

	//	_massIndepTuple = new TNtuple("_massIndepTuple","results from mass indep. fits", ("mass:eta:etaErr:f0:f0Err:f1:f1Err:f2:f2Err:eta2:eta2Err:etaReal:etaImg")   );
}

void JpsiGamEtaPiPiHist::plotDalitz(TH2F* theHisto,  EvtData* theData, double weight){
/*	Vector4<double>& V4_KsKl_HeliPsi = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_KsKl_HeliPsi]  ;
	Vector4<double>& V4_KpKm_HeliPsi = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_KpKm_HeliPsi] ;
	Vector4<double>& V4_gamma_HeliPsi = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_gamma_HeliPsi] ;

	double gphi1 = (V4_gamma_HeliPsi + V4_KsKl_HeliPsi ).M();
	double gphi2 = (V4_gamma_HeliPsi + V4_KpKm_HeliPsi ).M();
	theHisto->Fill( gphi1*gphi1, gphi2*gphi2   ,weight);*/
}

void JpsiGamEtaPiPiHist::plotEtaPipPim(TH1F* theHisto, EvtData* theData, double weight){
	Vector4<double>& v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPipPim_HeliPsi];
	theHisto->Fill( v4.M(), weight );
}

void JpsiGamEtaPiPiHist::plotEtaPi(TH1F* theHisto, EvtData* theData, double weight){
	Vector4<double>& v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi];
	theHisto->Fill( v4.M(), weight );
	v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi];
	theHisto->Fill( v4.M(), weight );
}
void JpsiGamEtaPiPiHist::plotPipPim(TH1F* theHisto, EvtData* theData, double weight){
	Vector4<double>& v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi];
	theHisto->Fill( v4.M(), weight );
}

void JpsiGamEtaPiPiHist::plotCostPhiEta(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtData* theData, double weight){
	Vector4<double>& v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliEtaPip];
	theCostHisto->Fill( v4.CosTheta(), weight );
	thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamEtaPiPiHist::plotCostPhiPip(TH1F* theCostHisto,  TH1F* thePhiHisto, EvtData* theData, double weight){
	Vector4<double>& v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPipPim];
	theCostHisto->Fill( v4.CosTheta(), weight );
	thePhiHisto->Fill( v4.Phi(), weight );
}

void JpsiGamEtaPiPiHist::plotCostGam(TH1F* theCostHisto, EvtData* theData, double weight){
	Vector4<double>& v4 = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_gamma_HeliPsi];
	theCostHisto->Fill( v4.CosTheta(), weight );
}
/*
void JpsiGamEtaPiPiHist::plotCostPhi_PhiPhiHeli(TH1F* theCostHisto, TH1F* thePhiHisto, const Vector4<double>& the4Vec, double weight){
	theCostHisto->Fill( the4Vec.CosTheta(), weight);
	thePhiHisto->Fill( the4Vec.Phi(), weight);
}
*/

void JpsiGamEtaPiPiHist::plotChi(TH1F* theChiHisto, EvtData* theData, double weight){
//   Vector4<double>& V4_KsKlKpKm_HeliPsi = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_KsKlKpKm_HeliPsi]  ;
//   Vector4<double>& V4_Ks_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Ks_HeliPsi] ;
//   Vector4<double>& V4_Kl_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Kl_HeliPsi] ;
//   Vector4<double>& V4_Kp_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Kp_HeliPsi] ;
//   Vector4<double>& V4_Km_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Km_HeliPsi] ;

//   double thePhiPhiDecayPlaneAngle = decayAngleChi( V4_KsKlKpKm_HeliPsi, V4_Kp_HeliPsi, V4_Km_HeliPsi, V4_Ks_HeliPsi, V4_Kl_HeliPsi   );


	Vector4<double>& V4_normPipPimDecHeliEtaPipPim = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_normPipPimDecHeliEtaPipPim];
	Vector4<double>& V4_normEtaPipDecHeliEtaPipPim = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_normEtaPipDecHeliEtaPipPim];

	double cosChi=(V4_normPipPimDecHeliEtaPipPim.Px()*V4_normEtaPipDecHeliEtaPipPim.Px()
				+V4_normPipPimDecHeliEtaPipPim.Py()*V4_normEtaPipDecHeliEtaPipPim.Py()
				+V4_normPipPimDecHeliEtaPipPim.Pz()*V4_normEtaPipDecHeliEtaPipPim.Pz())
				/ (V4_normPipPimDecHeliEtaPipPim.P()*V4_normEtaPipDecHeliEtaPipPim.P());

	double chi=acos(fabs(cosChi));
	theChiHisto->Fill(chi*180./TMath::Pi(),weight);
}

void  JpsiGamEtaPiPiHist::fillTuple( TNtuple* theTuple, EvtData* theData, double weight){
	Vector4<double>& V4_EtaPipPim_HeliPsi = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPipPim_HeliPsi]  ;
	Vector4<double>& V4_EtaPip_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi] ;
	Vector4<double>& V4_PipPim_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi] ;
	Vector4<double>& V4_gamma_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_gamma_HeliPsi] ;
  
	Vector4<double>& V4_EtaPip_HeliEtaPipPim= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliEtaPipPim] ;
	Vector4<double>& V4_PipPim_HeliKsKlKpKm= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliEtaPipPim] ;
	Vector4<double>& V4_Eta_HeliEtaPip= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliEtaPip] ;
	Vector4<double>& V4_Pip_HeliPipPim= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPipPim] ;
  
//	Vector4<double>& V4_Eta_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Eta_HeliPsi] ;
//	Vector4<double>& V4_Pip_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pip_HeliPsi] ;
//	Vector4<double>& V4_Pim_HeliPsi= theData->FourVecs[enumJpsiGamEtaPiPiData::V4_Pim_HeliPsi] ;
  
  
	Vector4<double>& V4_normPipPimDecHeliEtaPipPim = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_normPipPimDecHeliEtaPipPim];
	Vector4<double>& V4_normEtaPipDecHeliEtaPipPim = theData->FourVecs[enumJpsiGamEtaPiPiData::V4_normEtaPipDecHeliEtaPipPim];

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
			V4_PipPim_HeliKsKlKpKm.CosTheta(),
			V4_Eta_HeliEtaPip.CosTheta(),
			V4_Pip_HeliPipPim.CosTheta(),
			thePhiPhiDecayPlaneAngle,
//			testHeli,
			weight
			);
}


double JpsiGamEtaPiPiHist::decayAngleChi(const Vector4<double>& v4_p,const Vector4<double>& v4_d1,
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



