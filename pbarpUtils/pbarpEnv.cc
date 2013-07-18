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

// pbarpEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>

#include "pbarpUtils/pbarpEnv.hh"
#include "pbarpUtils/pbarpParser.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
//#include "PwaUtils/IsobarDecayList.hh"
#include "pbarpUtils/pbarpReaction.hh"
//#include "pbarpUtils/pbarpEventList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"

pbarpEnv* pbarpEnv::_instance=0;

pbarpEnv* pbarpEnv::instance()
{
  if (0==_instance) _instance = new pbarpEnv();
  return _instance;
}

pbarpEnv::pbarpEnv() :
  AbsEnv()
  ,_lmax(0)
  ,_pbarMomentum(0)
{
}
pbarpEnv::~pbarpEnv(){
}

void pbarpEnv::setup(pbarpParser* thePbarpParser){

  AbsEnv::setup(thePbarpParser);

  //Antiproton momentum
  _pbarMomentum = thePbarpParser->getpbarMomentum();

  double pMass=_particleTable->particle("proton")->mass();
  double antipMass=_particleTable->particle("antiproton")->mass();
  _initial4Vec = Vector4<double>(pMass+sqrt(antipMass*antipMass+_pbarMomentum*_pbarMomentum), 0., 0., _pbarMomentum);

  //Lmax
  _lmax=thePbarpParser->getLMax();

  // individual Lmax settings
  std::vector<std::string> theDropPbarpLForParticles = thePbarpParser->dropPbarpLForParticle();
  for(auto ldropIt = theDropPbarpLForParticles.begin(); ldropIt != theDropPbarpLForParticles.end(); ++ldropIt){
     std::string particle;
     short l;
     std::stringstream stringStr;
     stringStr << (*ldropIt);
     stringStr >> particle >> l;
     _dropPbarpLForParticleData[particle].push_back(l);
  }


  std::vector<std::string>::const_iterator itStr;


  //pbarp reaction
  _pbarpReaction=std::shared_ptr<pbarpReaction>(new pbarpReaction(_producedParticlePairs, _lmax));

  //fill prodDecayList
  if(thePbarpParser->productionFormalism()=="Cano"){
    std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs= _pbarpReaction->productionDecays();
    std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      _prodDecList->addDecay(*itDec);
    }
  }
  else if(thePbarpParser->productionFormalism()=="Tensor"){
      std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs= _pbarpReaction->productionTensorDecays();
      std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	 _prodDecList->addDecay(*itDec);
      }
  }
  else if(thePbarpParser->productionFormalism()=="Heli"){
     std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _pbarpReaction->productionHeliDecays();
     std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
     for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	_prodDecList->addDecay(*itDec);
     }
  }
  else{
    Alert <<"production formalism\t" << thePbarpParser->productionFormalism() << "\t is not supported!!!" << endmsg;
    exit(0);
  }

  //set prefactor for production and decay amplitudes
  std::map<std::string, double>::iterator strDoubleIt;
  for(strDoubleIt=_preFactorMap.begin(); strDoubleIt!=_preFactorMap.end(); ++strDoubleIt){
    std::string currentAmplitudeName=strDoubleIt->first;
    double currentPrefactor=strDoubleIt->second;
    
    std::shared_ptr<AbsDecay> currentDec=_prodDecList->decay(currentAmplitudeName);
    if(0!=currentDec){
      // currentDec->setPreFactor(currentPrefactor);
      // Info << "Set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
      currentDec->disableIsospin();
      Info << "Disable isospin coupling and set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
      
    }
    else{
      // look in decay amplitudes
      currentDec=_absDecList->decay(currentAmplitudeName);
      if(0!=currentDec){
	currentDec->setPreFactor(currentPrefactor);
	Info << "Set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg; 
      }
      else{
	Alert << "Amplitude with name\t" << currentAmplitudeName << "\tnot found!!!" << endmsg;
	exit(0);
      }
    }
  } 


  //set suffixes
  std::vector<std::string> suffixVec = thePbarpParser->replaceSuffixNames();
  std::map<std::string, std::string> decSuffixNames;

  for ( itStr = suffixVec.begin(); itStr != suffixVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string classStr;
    stringStr >> classStr;

    std::string suffixStr;
    stringStr >> suffixStr;
    decSuffixNames[classStr]=suffixStr;
  }
  
  //set suffixes for decay classes
  std::map<std::string, std::string>::iterator itMapStrStr;
  for (itMapStrStr=decSuffixNames.begin(); itMapStrStr!=decSuffixNames.end(); ++itMapStrStr){
    _absDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    _prodDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
    //    std::shared_ptr<IsobarDecay> theDec=_decList->decay(itMapStrStr->first);
  }

  //replace mass key
  std::vector<std::string> replMassKeyVec = thePbarpParser->replaceMassKey();
  std::map<std::string, std::string> decRepMassKeyNames;

  for ( itStr = replMassKeyVec.begin(); itStr != replMassKeyVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string oldStr;
    stringStr >> oldStr;

    std::string newStr;
    stringStr >> newStr;
    decRepMassKeyNames[oldStr]=newStr;
  }

  for (itMapStrStr=decRepMassKeyNames.begin(); itMapStrStr!=decRepMassKeyNames.end(); ++itMapStrStr){
    _absDecList->replaceMassKey(itMapStrStr->first, itMapStrStr->second);
  }


  //add dynamics 
  std::vector<std::shared_ptr<AbsDecay> > absDecList= _absDecList->getList();
  std::vector<std::string> decDynVec = thePbarpParser->decayDynamics();
  for ( itStr = decDynVec.begin(); itStr != decDynVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string particleStr;
    stringStr >> particleStr;

    std::string dynStr;
    stringStr >> dynStr;

    std::string tmpName;
    std::vector<std::string> additionalStringVec;
    while(stringStr >> tmpName){
      additionalStringVec.push_back(tmpName);
    }

    std::vector<std::shared_ptr<AbsDecay> >::iterator itDec;
    for (itDec=absDecList.begin(); itDec!=absDecList.end(); ++itDec){
      std::string theDecName=(*itDec)->name();
      std::string toFind=particleStr+"To";
      size_t found;
      found = theDecName.find(toFind);
      if (found!=string::npos && found==0){
	(*itDec)->enableDynamics(dynStr, additionalStringVec);
      }
    }
  }

  // spin density particles
  _spinDensity = thePbarpParser->spinDensityNames();

}


