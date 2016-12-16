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

// GammapChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf, Julian Pychy

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "gammapUtils/GammapChannelEnv.hh"
#include "ConfigParser/gammapParser.hh"
#include "gammapUtils/gammapReaction.hh"
#include "gammapUtils/gammapHist.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "ErrLogger/ErrLogger.hh"





GammapChannelEnv::GammapChannelEnv(gammapParser* theParser) : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_GAMMAP)
  ,_lmax(0)
  ,_theGamPParser(theParser)
{
}

void GammapChannelEnv::setup(ChannelID id){

   AbsChannelEnv::setup(id);

   //Lmax
   _lmax=_theGamPParser->getLMax();

   // individual Lmax settings
   std::vector<std::string> theDropGammapLForParticles = _theGamPParser->dropGammapLForParticle();
   for(auto ldropIt = theDropGammapLForParticles.begin(); ldropIt != theDropGammapLForParticles.end(); ++ldropIt){
      std::string particle;
      short l;
      std::stringstream stringStr;
      stringStr << (*ldropIt);
      stringStr >> particle >> l;
      _dropGammapLForParticleData[particle].push_back(l);
   }


   std::vector<std::string>::const_iterator itStr;


   //gammap reaction
   _gammapReaction=std::shared_ptr<gammapReaction>(new gammapReaction(_producedParticlePairs, id,_lmax));

   //fill prodDecayList
   std::vector<std::string> additionalStringVecDummy;
   std::string dynTypeDefault="WoDynamics";
   
   if(_theGamPParser->productionFormalism()=="Cano"){
      std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs= _gammapReaction->productionDecays();
      std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	if(_theGamPParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theGamPParser->qRProduction());
	else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
	_prodDecList->addDecay(*itDec);
      }
   }
   else if(_theGamPParser->productionFormalism()=="Tensor"){
      std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs= _gammapReaction->productionTensorDecays();
      std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	(*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
         _prodDecList->addDecay(*itDec);
      }
   }
   else if(_theGamPParser->productionFormalism()=="Heli"){
      std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _gammapReaction->productionHeliDecays();
      std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	if(_theGamPParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theGamPParser->qRProduction());
	else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
	_prodDecList->addDecay(*itDec);
      }
   }
   else{
      Alert <<"production formalism\t" << _theGamPParser->productionFormalism() << "\t is not supported!!!" << endmsg;
      exit(0);
   }


   //set prefactor for production and decay amplitudes
   std::map<std::string, double>::iterator strDoubleIt;
   for(strDoubleIt=_preFactorMap.begin(); strDoubleIt!=_preFactorMap.end(); ++strDoubleIt){
      std::string currentAmplitudeName=strDoubleIt->first;
      double currentPrefactor=strDoubleIt->second;

      std::shared_ptr<AbsDecay> currentDec=_prodDecList->decay(currentAmplitudeName);
      if(0!=currentDec){
	 currentDec->setPreFactor(currentPrefactor);
	 // InfoMsg << "Set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
	 currentDec->disableIsospin();
	 InfoMsg << "Disable isospin coupling and set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;

      }
      else{
	 // look in decay amplitudes
	 currentDec=_absDecList->decay(currentAmplitudeName);
	 if(0!=currentDec){
	    currentDec->setPreFactor(currentPrefactor);
	    InfoMsg << "Set prefactor " << currentPrefactor << " for amplitude " << currentAmplitudeName << endmsg;
	 }
	 else{
	    Alert << "Amplitude with name\t" << currentAmplitudeName << "\tnot found!!!" << endmsg;
	    exit(0);
	 }
      }
   }



   //set suffixes
   std::vector<std::string> suffixVec = _theGamPParser->replaceSuffixNames();
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
   }


   //replace mass key
   std::vector<std::string> replMassKeyVec = _theGamPParser->replaceMassKey();
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
   std::vector<std::string> decDynVec = _theGamPParser->decayDynamics();
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
     (*itProdDecList)->setDecayLevelTree(AbsDecay::decayLevel::isProdAmp);    
   }
   
   // spin density particles
   _spinDensity = _theGamPParser->spinDensityNames();

}



std::shared_ptr<AbsHist> GammapChannelEnv::CreateHistInstance(std::string additionalSuffix){

  return std::shared_ptr<AbsHist>(new gammapHist(additionalSuffix));
}
