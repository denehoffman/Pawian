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
#include "PwaUtils/ProdChannelInfo.hh"
#include "ErrLogger/ErrLogger.hh"





EpemChannelEnv::EpemChannelEnv(epemParser* theParser) : 
  AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_EPEM)
  ,_theEpEmParser(theParser)
  ,_cmsMass(1.)
{
}

void EpemChannelEnv::setup(ChannelID id){

   AbsChannelEnv::setup(id);


   _cmsMass=_theEpEmParser->cmsMass();

  // has to be set via parser !!!!
  double totalyMom=0.04;
  _initial4Vec = Vector4<double>( sqrt(_cmsMass*_cmsMass+totalyMom*totalyMom), 0., totalyMom, 0.);

  std::vector<std::string>::const_iterator itStr;


  //epem reaction
  _epemReaction=std::shared_ptr<epemReaction>(new epemReaction(_prodChannelInfoList, id));
  std::vector<std::string> additionalStringVecDummy;
  std::string dynTypeDefault="WoDynamics";

  //  std::vector< std::shared_ptr<AbsDecay> > prodDecs;
  if (_theEpEmParser->productionFormalism()=="Heli" || _theEpEmParser->productionFormalism()=="HeliMultipole"){
    std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs;
    if (_theEpEmParser->productionFormalism()=="Heli") prodDecs = _epemReaction->productionHeliDecays();
    else prodDecs = _epemReaction->productionHeliMultipoleDecays();
    std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
  else if (_theEpEmParser->productionFormalism()=="Tensor"){
    std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs;
    if(_theEpEmParser->productionTensorRadType() == "Zou"){
      Info << "use productionTensorRadType " << _theEpEmParser->productionTensorRadType() << endmsg;
      prodDecs=_epemReaction->productionTensorZouDecays();
    } 
    else if(_theEpEmParser->productionTensorRadType() == "Default"){
      Info << "use productionTensorRadType " << _theEpEmParser->productionTensorRadType() << endmsg;
      prodDecs=_epemReaction->productionTensorDecays(); //default
    } 
    else{
      Alert <<"productionTensorRadType with the name " << _theEpEmParser->productionTensorRadType() << " doesn't exist!!!" << endmsg;
      exit(0);
    }


    std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
    for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
      if((*itDec)->prodChannelInfo()->withProdBarrier()) (*itDec)->enableProdBarrier();
      else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      _prodDecList->addDecay(*itDec);
    }
  }
   else{
      Alert <<"production formalism\t" << _theEpEmParser->productionFormalism() << "\t is not supported!!!" << endmsg;
      exit(0);
   }

  // std::vector< std::shared_ptr<AbsDecay> >::iterator itDec;
  // for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
  //   _prodDecList->addDecay(*itDec);
  // }

  //set prefactor for production and decay amplitudes
  std::map<std::string, double>::iterator strDoubleIt;
  for(strDoubleIt=_preFactorMap.begin(); strDoubleIt!=_preFactorMap.end(); ++strDoubleIt){
      std::string currentAmplitudeName=strDoubleIt->first;
      double currentPrefactor=strDoubleIt->second;

      std::shared_ptr<AbsDecay> currentDec=_prodDecList->decay(currentAmplitudeName);
      if(0!=currentDec){
         currentDec->setPreFactor(currentPrefactor);
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
  std::vector<std::string> suffixVec = _theEpEmParser->replaceSuffixNames();
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
  std::vector<std::string> replMassKeyVec = _theEpEmParser->replaceMassKey();
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
  std::vector<std::string> decDynVec = _theEpEmParser->decayDynamics();
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
     currentProdAmp->setDecayLevelTree(AbsDecay::decayLevel::isProdAmp, currentProdAmp, currentProdAmp);    
   }
}



std::shared_ptr<AbsHist> EpemChannelEnv::CreateHistInstance(std::string additionalSuffix){

  return std::shared_ptr<AbsHist>(new epemHist(additionalSuffix));
}
