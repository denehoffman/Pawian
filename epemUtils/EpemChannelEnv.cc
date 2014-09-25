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

// EpemChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "epemUtils/EpemChannelEnv.hh"
#include "ConfigParser/epemParser.hh"
#include "epemUtils/epemReaction.hh"
#include "epemUtils/epemHist.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "ErrLogger/ErrLogger.hh"





EpemChannelEnv::EpemChannelEnv(epemParser* theParser) : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_EPEM)
  ,_theParser(theParser)
{
}

void EpemChannelEnv::setup(ChannelID id){

   AbsChannelEnv::setup(id);


   _cmsMass=_theParser->cmsMass();

  // has to be set via parser !!!!
  double totalyMom=0.04;
  _initial4Vec = Vector4<double>( sqrt(_cmsMass*_cmsMass+totalyMom*totalyMom), 0., totalyMom, 0.);

  std::vector<std::string>::const_iterator itStr;


  //epem reaction
  _epemReaction=std::shared_ptr<epemReaction>(new epemReaction(_producedParticlePairs, id));
  std::vector<std::string> additionalStringVecDummy;

  std::vector< std::shared_ptr<AbsDecay> > prodDecs;
  if (_theParser->productionFormalism()=="Heli"){
    std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _epemReaction->productionHeliDecays();
    std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if(_theParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theParser->qRProduction());
      _prodDecList->addDecay(*itDec);
    }
  }
  else if (_theParser->productionFormalism()=="Tensor"){
    std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs= _epemReaction->productionTensorDecays();
    std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if(_theParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theParser->qRProduction());
      _prodDecList->addDecay(*itDec);
    }
  }
   else{
      Alert <<"production formalism\t" << _theParser->productionFormalism() << "\t is not supported!!!" << endmsg;
      exit(0);
   }

  std::vector< std::shared_ptr<AbsDecay> >::iterator itDec;
  for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
    _prodDecList->addDecay(*itDec);
  }

  //set suffixes
  std::vector<std::string> suffixVec = _theParser->replaceSuffixNames();
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
  std::vector<std::string> replMassKeyVec = _theParser->replaceMassKey();
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
  std::vector<std::string> decDynVec = _theParser->decayDynamics();
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
}



void EpemChannelEnv::CreateHistInstance(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams, std::string additionalSuffix){

   epemHist(theLh, theFitParams, additionalSuffix);
}
