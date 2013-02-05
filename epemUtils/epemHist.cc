// epemHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "epemUtils/epemHist.hh"
#include "epemUtils/epemEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/EvtDataBaseList.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

epemHist::epemHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
  AbsHist() 
{
  initRootStuff();
  fillIt(theLh, theFitParams);
}

epemHist::~epemHist(){
  _theTFile->Write();
  _theTFile->Close();
}

void epemHist::initRootStuff(){
  std::ostringstream rootFileName;
  rootFileName << "./pawianHists" << epemEnv::instance()->outputFileNameSuffix() << ".root";
  _theTFile=new TFile(rootFileName.str().c_str(),"recreate");

  std::vector<std::vector<std::string> > histMassNameVec=epemEnv::instance()->histMassSystems();
  std::vector<std::vector<std::string> >::iterator itVecStr;
  for(itVecStr=histMassNameVec.begin(); itVecStr!=histMassNameVec.end(); ++itVecStr){
    boost::shared_ptr<massHistData> tmpMassHistData(new massHistData(*itVecStr));
    std::string tmpBaseName=tmpMassHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histName="data"+tmpBaseName;
    std::string histDescription = "M("+tmpMassHistData->_name+") (data)";

    double psiMass = epemEnv::instance()->particleTable()->particle("Jpsi")->mass();
    double massMin = 0;
    double massMax = psiMass;

    std::vector<std::string> fspNames=tmpMassHistData->_fspNames;
    std::vector<Particle*> allFsp =  epemEnv::instance()->finalStateParticles();
    std::vector<Particle*>::iterator itAllFsp;

    for(itAllFsp = allFsp.begin(); itAllFsp != allFsp.end(); ++itAllFsp){
       bool isObserver = true;
       std::vector<std::string>::iterator itStr2;
       for(itStr2=fspNames.begin(); itStr2!=fspNames.end(); ++itStr2){
	  if(*itStr2 == (*itAllFsp)->name())
	     isObserver = false;
       }
       if(isObserver)
	  massMax -= (*itAllFsp)->mass();
       else
	  massMin += (*itAllFsp)->mass();
    }

    massMax += (massMax - massMin) * 0.02;
    massMin -= (massMax - massMin) * 0.02;

    TH1F* currentMassDataHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
    currentMassDataHist->Sumw2();
    _massDataHistMap[tmpMassHistData]=currentMassDataHist;

    histName="MC"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (MC)";
    TH1F* currentMassMcHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
    currentMassMcHist->Sumw2();
    _massMcHistMap[tmpMassHistData]=currentMassMcHist;

    histName="Fit"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (fit)";
    TH1F* currentMassFitHist=new TH1F(histName.c_str(), histDescription.c_str(), 100., massMin, massMax);
    currentMassFitHist->Sumw2();
    _massFitHistMap[tmpMassHistData]=currentMassFitHist;
  }

  std::vector<boost::shared_ptr<angleHistData> > angleHistDataVec=epemEnv::instance()->angleHistDataVec();

  std::vector<boost::shared_ptr<angleHistData> >::iterator itAngleVec;
  for (itAngleVec=angleHistDataVec.begin(); itAngleVec!=angleHistDataVec.end(); ++itAngleVec){
    std::string tmpBaseName= (*itAngleVec)->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histThetaName="dataTheta"+tmpBaseName;
    std::string histPhiName="dataPhi"+tmpBaseName;
    std::string histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (data)";
    std::string histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (data)";

    TH1F* currentThetaAngleDataHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleDataHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleDataHist->Sumw2();
    currentPhiAngleDataHist->Sumw2();
    _angleDataHistMap[*itAngleVec]=std::pair<TH1F*, TH1F*>(currentThetaAngleDataHist, currentPhiAngleDataHist);

    histThetaName="MCTheta"+tmpBaseName;
    histPhiName="MCPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (MC)";
    histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (MC)";

    TH1F* currentThetaAngleMcHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleMcHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleMcHist->Sumw2();
    currentPhiAngleMcHist->Sumw2();
    _angleMcHistMap[*itAngleVec]=std::pair<TH1F*, TH1F*>(currentThetaAngleMcHist, currentPhiAngleMcHist);

    histThetaName="FitTheta"+tmpBaseName;
    histPhiName="FitPhi"+tmpBaseName;
    histThetaDescription = "cos#Theta(" +(*itAngleVec)->_name + ") (fit)";
    histPhiDescription = "#phi(" +(*itAngleVec)->_name + ") (fit)";

    TH1F* currentThetaAngleFitHist=new TH1F(histThetaName.c_str(), histThetaDescription.c_str(), 100., -1., 1.);
    TH1F* currentPhiAngleFitHist=new TH1F(histPhiName.c_str(), histPhiDescription.c_str(), 100., -3.14159, 3.14159);
    currentThetaAngleFitHist->Sumw2();
    currentPhiAngleFitHist->Sumw2();
    _angleFitHistMap[*itAngleVec]=std::pair<TH1F*, TH1F*>(currentThetaAngleFitHist, currentPhiAngleFitHist);
  } 

 std::map<boost::shared_ptr<angleHistData>, TH1F*, pawian::Collection::SharedPtrLess > _angleDataHistMap;
}

