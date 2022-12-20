//************************************************************************//
//									  //
//  Copyright 20177Bertram Kopf (bertram@ep1.rub.de)			  //
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

// PiPiScatteringChannelEnv class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include "Utils/PawianIOUtils.hh"
#include "Particle/ParticleTable.hh"
#include "pipiScatteringUtils/PiPiScatteringChannelEnv.hh"
#include "ConfigParser/pipiScatteringParser.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/PiPiScatteringDecay.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "PwaUtils/FsParticleProjections.hh"
#include "PwaUtils/RootPiPiScatteringHist.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Event/MassRangeCut.hh"

PiPiScatteringChannelEnv::PiPiScatteringChannelEnv(pipiScatteringParser* theParser) 
  : AbsChannelEnv(theParser, AbsChannelEnv::CHANNEL_PIPISCATTERING)
{
}



void PiPiScatteringChannelEnv::setupChannel(ChannelID id){

  PiPiScatteringChannelEnv::setupGlobal(id);

  //set prefactor for production and decay amplitudes
   //   AbsChannelEnv::setPrefactors();
  
   //replace suffixes for fit parameter
  AbsChannelEnv::replaceParameterSuffixes();
 
  //replace mass key
  AbsChannelEnv::replaceMassKeys();

  //add dynamics
  //AbsChannelEnv::addDynamics();
  addDynamics();
}

void PiPiScatteringChannelEnv::setupGlobal(ChannelID id){
   if(_alreadySetUp){
      Alert << "PiPiScatteringChannelEnv already set up!" << endmsg;
      exit(1);
   }
   _alreadySetUp = true;

   // Set channel id
   _channelID = id;


   // cloned particles
   const std::vector<std::string> cloneParticle=_theParser->cloneParticle();
   std::vector<std::string>::const_iterator itcP;

   for ( itcP = cloneParticle.begin(); itcP != cloneParticle.end(); ++itcP){
      std::istringstream particles(*itcP);
      std::string particleOld;
      std::string particleNew;
      particles >> particleOld >> particleNew;

      GlobalEnv::instance()->particleTable()->clone(particleNew, particleOld);
   }

   //final state particles
   const std::vector<std::string> finalStateParticleStr=_theParser->finalStateParticles();
   
   std::vector<std::string>::const_iterator itStr;
   for ( itStr = finalStateParticleStr.begin(); itStr != finalStateParticleStr.end(); ++itStr){
     Particle* currentParticle = GlobalEnv::instance()->particleTable()->particle(*itStr);
     _finalStateParticles.push_back(currentParticle);
   }
   
   _noFinalStateParticles= (int) _finalStateParticles.size();
   if (_noFinalStateParticles !=2){
     Alert << "number of final state particles must be exactly 2!!!" 
	   << "\nhere: " << _noFinalStateParticles << " final state particles" << endmsg;
     exit(1);
   }

   if( (_theParser->productionFormalism()!="PVecIntensity")
          && (_theParser->productionFormalism()!="PhaseDiff") ){
     InfoMsg << "ProForm: " << _theParser->productionFormalism() << endmsg;
     _fsParticleProjections = std::shared_ptr<FsParticleProjections>(new FsParticleProjections(_finalStateParticles));
   }
   
   //decays
   std::vector<std::string> decaySystem= _theParser->decaySystem();
   if (decaySystem.size() !=1){
     Alert << "not more than 1 decay is required!!!" 
	   << "here: " << decaySystem.size() << " decays have been set up" << endmsg;
     exit(1);
   }

   std::stringstream stringStr;
   stringStr << decaySystem.at(0);

   Particle* motherParticle =0;
   std::vector<Particle*> daughterParticles;
   std::string tmpName;
   bool isDecParticle=false;

   while(stringStr >> tmpName){
     if(tmpName=="To") {
       isDecParticle=true;
       continue;
     }
     if (GlobalEnv::instance()->particleTable()->particle(tmpName) == NULL){
       Alert << "Particle " << tmpName << " does not exist in ParticleTable. Please clone an existing particle or add it to the table.";
       exit(1);
     }
     if(isDecParticle){
       daughterParticles.push_back(GlobalEnv::instance()->particleTable()->particle(tmpName));
     }
     else{
       motherParticle = GlobalEnv::instance()->particleTable()->particle(tmpName);
     }
   }   

   if(daughterParticles.size()!=2){
     Alert << "it is required to define exactly 2 daughter particles!!!"
	   <<"\nths number of daughter particles here are " << daughterParticles.size() 
	   << endmsg;
     exit(1);
   }

   std::shared_ptr<AbsDecay> tmpDec= std::shared_ptr<AbsDecay>(new PiPiScatteringDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));;

   _absDecList->addDecay(tmpDec);
    InfoMsg << "added decay " << tmpDec->name() << " to decay list" << endmsg;


  //create global map for fit parameter suffixes to be replaced
  std::vector<std::string> suffixVec = _theParser->replaceSuffixNames();
  for ( itStr = suffixVec.begin(); itStr != suffixVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string classStr;
    stringStr >> classStr;

    std::string suffixStr;
    stringStr >> suffixStr;

    GlobalEnv::instance()->addIntoToBeReplacedSuffixMap(classStr, suffixStr);
  }

  _channelSubTypeName=motherParticle->name()+_theParser->productionFormalism();


  //mass range cuts
  std::vector<std::string> massRangeCuts=_theParser->massRangeCuts();
  if(massRangeCuts.size()>0) _useMassRange=true;
  for ( itStr = massRangeCuts.begin(); itStr != massRangeCuts.end(); ++itStr){
    std::string currentString=*itStr;
    std::shared_ptr<MassRangeCut> currentMassRangeCut(new MassRangeCut( currentString, _finalStateParticles));
    _massRangeCuts.push_back(currentMassRangeCut);
  }
}

void PiPiScatteringChannelEnv::addDynamics(){
  std::vector<std::shared_ptr<AbsDecay> > absDecList= _absDecList->getList();
  std::vector<std::string> decDynVec = _theParser->decayDynamics();
  
  if (decDynVec.size() !=1){
    Alert << "not more than 1 description for the dynamics is required!!!" 
	  << "here: " << decDynVec.size() << " dyn types have been set up" << endmsg;
    exit(1);
  }
  
  std::string dynCfgLine=decDynVec.at(0);
  
  std::stringstream stringStr;
  stringStr << dynCfgLine;

  std::string particleStr;
  stringStr >> particleStr;

  std::string dynStr;
  stringStr >> dynStr;

  std::string tmpName;
  std::vector<std::string> additionalStringVec;
  while(stringStr >> tmpName){
    additionalStringVec.push_back(tmpName);
  }

  InfoMsg << "dynStr: " << dynStr << endmsg;

  if(dynStr=="KMatrix" || dynStr=="TMatrix" || dynStr=="FVectorIntensity"|| dynStr=="TMatrixCompare" || dynStr=="FVectorCompare") _pathKMatrixParserFile=additionalStringVec[0];
  else{
    Alert << "dyn type " << dynStr << " is not supported for pipi scattering reactions!!!" << endmsg;
    exit(0);
  }

  _pathKMatrixParserFile=PawianIOUtils::getFileName(GlobalEnv::instance()->KMatrixStorePath(), _pathKMatrixParserFile);

  std::shared_ptr<AbsDecay> theDecay=absDecList.at(0);

  theDecay->enableDynamics(dynStr, additionalStringVec);

}

std::shared_ptr<AbsHist> PiPiScatteringChannelEnv::CreateHistInstance(std::string additionalSuffix, bool withTruth){
  return std::shared_ptr<AbsHist>(new RootPiPiScatteringHist(additionalSuffix, withTruth));
}



