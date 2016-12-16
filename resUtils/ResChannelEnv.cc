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

// ResChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "Particle/ParticleTable.hh"
#include "resUtils/ResChannelEnv.hh"
#include "ConfigParser/resParser.hh"
#include "resUtils/resHist.hh"
#include "resUtils/resReaction.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

ResChannelEnv::ResChannelEnv(resParser* theParser) 
  : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_RES)
  ,_theResParser(theParser)
  ,_polarizedMother(false)
{
}



void ResChannelEnv::setup(ChannelID id){

   AbsChannelEnv::setup(id);

   //  // has to be set via parser !!!!
  // double totalyMom=0.04;
  // _initial4Vec = Vector4<double>( sqrt(_cmsMass*_cmsMass+totalyMom*totalyMom), 0., totalyMom, 0.);

  _motherParticle = GlobalEnv::instance()->particleTable()->particle(_theResParser->motherResName());
  if(0==_motherParticle){
    Alert << "mother particle with name\t" << _theResParser->motherResName() << "\tdoesn't exist" << endmsg;
    exit(0);
  }

  _polarizedMother=_theResParser->polarizedMother();
  if(_polarizedMother) {
    InfoMsg << "polarization of the mother particle will be taken into accout!!!" << endmsg;
  }

  //  double theMotherMass=_motherParticle->mass();
  _initial4Vec = Vector4<double>(_motherParticle->mass(), 0. ,0. ,0.);
  std::vector<std::string>::const_iterator itStr;


  //epem reaction
  _resReaction=std::shared_ptr<resReaction>(new resReaction(_motherParticle, _prodChannelInfoList, id));

  //fill prodDecayList
  std::vector<std::string> additionalStringVecDummy;
  std::string dynTypeDefault="WoDynamics";

  if(_theResParser->productionFormalism()=="Heli"){
    std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _resReaction->productionHeliDecays();
    std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      (*itDec)->disableIsospin();
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else if(_theResParser->productionFormalism()=="Cano"){
    std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs= _resReaction->productionCanoDecays();
    std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      (*itDec)->disableIsospin();
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else{
    Alert << "production formalism: " << _theResParser->productionFormalism() << " is not supported for res reactions" << endmsg;
    exit(0); 
  }

  //set suffixes
  std::vector<std::string> suffixVec = _theResParser->replaceSuffixNames();
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
  std::vector<std::string> replMassKeyVec = _theResParser->replaceMassKey();
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
  std::vector<std::string> decDynVec = _theResParser->decayDynamics();
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

   //set decay levels
   std::vector<std::shared_ptr<AbsDecay> > prodDecList= _prodDecList->getList();
   std::vector<std::shared_ptr<AbsDecay> >::iterator itProdDecList;
   for (itProdDecList=prodDecList.begin(); itProdDecList!=prodDecList.end(); ++itProdDecList){
     std::shared_ptr<AbsDecay> currentProdAmp= (*itProdDecList);
     (*itProdDecList)->setDecayLevelTree(AbsDecay::decayLevel::isProdAmp, currentProdAmp, currentProdAmp);    
   }
}



std::shared_ptr<AbsHist> ResChannelEnv::CreateHistInstance(std::string additionalSuffix){

  return std::shared_ptr<AbsHist>(new resHist(additionalSuffix));
}
