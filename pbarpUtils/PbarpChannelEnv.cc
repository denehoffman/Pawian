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

// PbarpChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "pbarpUtils/PbarpChannelEnv.hh"
#include "ConfigParser/pbarpParser.hh"
#include "pbarpUtils/pbarpReaction.hh"
#include "pbarpUtils/pbarpHist.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "ErrLogger/ErrLogger.hh"





PbarpChannelEnv::PbarpChannelEnv(pbarpParser* theParser) : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_PBARP)
  ,_lmax(0)
  ,_pbarMomentum(0)
  ,_theParser(theParser)
{
}

void PbarpChannelEnv::setup(ChannelID id){

   AbsChannelEnv::setup(id);

   //Antiproton momentum
   _pbarMomentum = _theParser->getpbarMomentum();

   double pMass=GlobalEnv::instance()->particleTable()->particle("proton")->mass();
   double antipMass=GlobalEnv::instance()->particleTable()->particle("antiproton")->mass();
   _initial4Vec = Vector4<double>(pMass+sqrt(antipMass*antipMass+_pbarMomentum*_pbarMomentum), 0., 0., _pbarMomentum);

   //Lmax
   _lmax=_theParser->getLMax();

   // individual Lmax settings
   std::vector<std::string> theDropPbarpLForParticles = _theParser->dropPbarpLForParticle();
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
   _pbarpReaction=std::shared_ptr<pbarpReaction>(new pbarpReaction(_producedParticlePairs, id,_lmax));

   //preparations for prod key replacements
   std::vector<std::string> replProdKeyVec = _theParser->replaceProdKey();
   std::map<std::string, std::string> repProdKeyNames;

   for ( itStr = replProdKeyVec.begin(); itStr != replProdKeyVec.end(); ++itStr){
      std::stringstream stringStr;
      stringStr << (*itStr);
      std::string oldStr;
      stringStr >> oldStr;

      std::string newStr;
      stringStr >> newStr;
      repProdKeyNames[oldStr]=newStr;
   }


   // //preparations for mass key replacements
   // std::vector<std::string> replMassKeyVec = _theParser->replaceMassKey();
   // std::map<std::string, std::string> decRepMassKeyNames;

   // for ( itStr = replMassKeyVec.begin(); itStr != replMassKeyVec.end(); ++itStr){
   //    std::stringstream stringStr;
   //    stringStr << (*itStr);
   //    std::string oldStr;
   //    stringStr >> oldStr;

   //    std::string newStr;
   //    stringStr >> newStr;
   //    decRepMassKeyNames[oldStr]=newStr;
   // }

   //fill prodDecayList
   std::vector<std::string> additionalStringVecDummy;
   std::string dynTypeDefault="WoDynamics";

   additionalStringVecDummy.push_back("0.197");
   if(_theParser->productionFormalism()=="Cano"){
      std::vector< std::shared_ptr<IsobarLSDecay> > prodDecs= _pbarpReaction->productionDecays();
      std::vector< std::shared_ptr<IsobarLSDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	//	if(_theParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theParser->qRProduction());
	//	else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
	 _prodDecList->addDecay(*itDec);
      }
   }
   else if(_theParser->productionFormalism()=="Tensor"){
      std::vector< std::shared_ptr<IsobarTensorDecay> > prodDecs= _pbarpReaction->productionTensorDecays();
      std::vector< std::shared_ptr<IsobarTensorDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
         _prodDecList->addDecay(*itDec);
	 //	 (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
      }
   }
   else if(_theParser->productionFormalism()=="Heli"){
      std::vector< std::shared_ptr<IsobarHeliDecay> > prodDecs= _pbarpReaction->productionHeliDecays();
      std::vector< std::shared_ptr<IsobarHeliDecay> >::iterator itDec;
      for (itDec=prodDecs.begin(); itDec!=prodDecs.end(); ++itDec){
	//	if(_theParser->useProductionBarrier()) (*itDec)->enableProdBarrier(_theParser->qRProduction());
	//	else (*itDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
	 _prodDecList->addDecay(*itDec);
      }
   }
   else{
      Alert <<"production formalism\t" << _theParser->productionFormalism() << "\t is not supported!!!" << endmsg;
      exit(0);
   }

   //enable dynamics and replace mass keys
   std::map<std::string, std::string>::iterator itMapStrStr;
   for (itMapStrStr=repProdKeyNames.begin(); itMapStrStr!=repProdKeyNames.end(); ++itMapStrStr){
      _prodDecList->replaceProdKey(itMapStrStr->first, itMapStrStr->second);
   }

   std::vector< std::shared_ptr<AbsDecay> > theProdDecs=_prodDecList->getList();
   std::vector< std::shared_ptr<AbsDecay> >::iterator itAbsDec;
   for (itAbsDec=theProdDecs.begin(); itAbsDec!=theProdDecs.end(); ++itAbsDec){
     if(_theParser->useProductionBarrier()) (*itAbsDec)->enableProdBarrier(_theParser->qRProduction());
     else (*itAbsDec)->enableDynamics(dynTypeDefault, additionalStringVecDummy);
   }


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
   //   std::map<std::string, std::string>::iterator itMapStrStr;
   for (itMapStrStr=decSuffixNames.begin(); itMapStrStr!=decSuffixNames.end(); ++itMapStrStr){
      _absDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
      _prodDecList->replaceSuffix(itMapStrStr->first, itMapStrStr->second);
      //    std::shared_ptr<IsobarDecay> theDec=_decList->decay(itMapStrStr->first);
   }


   //replace mass key for decays
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


   //set decay levels
   std::vector<std::shared_ptr<AbsDecay> > prodDecList= _prodDecList->getList();
   std::vector<std::shared_ptr<AbsDecay> >::iterator itProdDecList;
   for (itProdDecList=prodDecList.begin(); itProdDecList!=prodDecList.end(); ++itProdDecList){
     (*itProdDecList)->setDecayLevelTree(AbsDecay::decayLevel::isProdAmp);    
   } 

   // spin density particles
   _spinDensity = _theParser->spinDensityNames();

}



void PbarpChannelEnv::CreateHistInstance(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams, std::string additionalSuffix){

   pbarpHist(theLh, theFitParams, additionalSuffix);
}
