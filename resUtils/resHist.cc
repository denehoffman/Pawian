//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

// resHist class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "resUtils/resHist.hh"
#include "resUtils/resEnv.hh"
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

resHist::resHist(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams) :
  AbsHist(resEnv::instance()) 
{
  initRootStuff();
  fillIt(theLh, theFitParams);
}

resHist::~resHist(){
}

void resHist::initRootStuff(){

  std::vector<std::vector<std::string> > histMassNameVec=resEnv::instance()->histMassSystems();
  std::vector<std::vector<std::string> >::iterator itVecStr;
  for(itVecStr=histMassNameVec.begin(); itVecStr!=histMassNameVec.end(); ++itVecStr){
    std::shared_ptr<massHistData> tmpMassHistData(new massHistData(*itVecStr));
    std::string tmpBaseName=tmpMassHistData->_name;
    boost::replace_all(tmpBaseName,"+","p");
    boost::replace_all(tmpBaseName,"-","m");
    std::string histName="Data"+tmpBaseName;
    std::string histDescription = "M("+tmpMassHistData->_name+") (data)";

    //    double psiMass = epemEnv::instance()->particleTable()->particle("Jpsi")->mass();
    double massMin = 0;
    //    double massMax = psiMass;
    double massMax = resEnv::instance()->motherParticle()->mass();
    //    DebugMsg << "epemEnv::instance()->cmsMass():\t" << epemEnv::instance()->cmsMass() << endmsg;

    std::vector<std::string> fspNames=tmpMassHistData->_fspNames;
    std::vector<Particle*> allFsp =  resEnv::instance()->finalStateParticles();
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

    TH1F* currentMassDataHist=new TH1F(histName.c_str(), histDescription.c_str(), 100, massMin, massMax);
    currentMassDataHist->Sumw2();
    _massDataHistMap[tmpMassHistData]=currentMassDataHist;

    DebugMsg << histName.c_str() << "\t" << massMin << "\t" << massMax << endmsg;

    histName="MC"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (MC)";
    TH1F* currentMassMcHist=new TH1F(histName.c_str(), histDescription.c_str(), 100, massMin, massMax);
    currentMassMcHist->Sumw2();
    _massMcHistMap[tmpMassHistData]=currentMassMcHist;

    DebugMsg << histName.c_str() << "\t" << massMin << "\t" << massMax << endmsg;
    
    histName="Fit"+tmpBaseName;
    histDescription = "M("+tmpMassHistData->_name+") (fit)";
    TH1F* currentMassFitHist=new TH1F(histName.c_str(), histDescription.c_str(), 100, massMin, massMax);
    currentMassFitHist->Sumw2();
    _massFitHistMap[tmpMassHistData]=currentMassFitHist;

    DebugMsg << histName.c_str() << "\t" << massMin << "\t" << massMax << endmsg;
  }

}

