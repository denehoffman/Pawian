// pbarpHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "pbarpUtils/pbarpHist.hh"
#include "pbarpUtils/pbarpEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
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
//#include "TMath.h"

pbarpHist::pbarpHist(boost::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
  AbsHist(pbarpEnv::instance()) 
{
  initRootStuff();
  fillIt(theLh, theFitParams);
}

pbarpHist::~pbarpHist(){
  // _theTFile->Write();
  // _theTFile->Close();
}

void pbarpHist::initRootStuff(){


  std::vector<std::vector<std::string> > histMassNameVec=pbarpEnv::instance()->histMassSystems();
  std::vector<std::vector<std::string> >::iterator itVecStr;
  for(itVecStr=histMassNameVec.begin(); itVecStr!=histMassNameVec.end(); ++itVecStr){
    boost::shared_ptr<massHistData> tmpMassHistData(new massHistData(*itVecStr));
    std::string tmpBaseName=tmpMassHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histName="Data"+tmpBaseName;
    std::string histDescription = "M("+tmpMassHistData->_name+") (data)";

    double pMass = pbarpEnv::instance()->particleTable()->particle("proton")->mass();
    double pbarMom = pbarpEnv::instance()->pbarMomentum();
    double massMin = 0;
    double massMax = sqrt(pow(sqrt(pMass*pMass + pbarMom*pbarMom) + pMass, 2) - pbarMom*pbarMom);

    std::vector<std::string> fspNames=tmpMassHistData->_fspNames;
    std::vector<Particle*> allFsp =  pbarpEnv::instance()->finalStateParticles();
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

}


