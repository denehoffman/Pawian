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

// AbsChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy


#include "AbsChannelEnv.hh"
#include "ConfigParser/ParserBase.hh"
#include "GlobalEnv.hh"
#include "AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "PwaUtils/IsobarTensorPsiToGamXDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "PwaUtils/AbsLh.hh"
#include "Particle/Particle.hh"
#include "Particle/PdtParser.hh" 
#include "Event/MassRangeCut.hh"
#include "ErrLogger/ErrLogger.hh"



short AbsChannelEnv::CHANNEL_PBARP = 1;
short AbsChannelEnv::CHANNEL_EPEM = 2;
short AbsChannelEnv::CHANNEL_RES = 3;
short AbsChannelEnv::CHANNEL_GAMMAP = 4;

AbsChannelEnv::AbsChannelEnv(ParserBase* theParser, short channelType) :
  _channelType(channelType)
  ,_alreadySetUp(false)
  ,_noFinalStateParticles(0)
  ,_absDecList(new AbsDecayList())
  ,_prodDecList(new AbsDecayList())
  ,_useMassRange(false)
  //  ,_massMin(0.)
  //  ,_massMax(100.)
  ,_theParser(theParser)
{
   _theLh.reset();
}



std::shared_ptr<AbsLh> AbsChannelEnv::Lh(){
   if(_theLh == nullptr){
      Alert << "Accessing Lh NULL pointer." << endmsg;
   }
   return _theLh;
}



void AbsChannelEnv::setup(ChannelID id){
   if(_alreadySetUp){
      Alert << "PbarpChannelEnv already set up!" << endmsg;
      exit(1);
   }
   _alreadySetUp = true;

   // Set channel id
   _channelID = id;

   // Event weights for data?
   _useEvtWeight = _theParser->useEvtWeight();

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

  //decays

  std::vector<std::string> decaySystem= _theParser->decaySystem();
  for ( itStr = decaySystem.begin(); itStr != decaySystem.end(); ++itStr){

    Particle* motherParticle =0;
    std::vector<Particle*> daughterParticles;

    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;

    bool useIsospin=true;
    bool isDecParticle=false;
    bool firstArgument=true;
    bool secondArgument=false;
    std::string usedSystem("");

    while(stringStr >> tmpName){
      if(firstArgument){
	usedSystem=tmpName;
	firstArgument=false;
	secondArgument=true;
	continue;
      }
      if(secondArgument){
	if(tmpName=="noIso"){
	 useIsospin=false;
         continue;
        }
	secondArgument=false;
      }
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
    std::shared_ptr<AbsDecay> tmpDec;
    if(daughterParticles.size()==2){
      if (usedSystem=="Heli") tmpDec= std::shared_ptr<AbsDecay>(new IsobarHeliDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
      else if (usedSystem=="Cano")  tmpDec= std::shared_ptr<AbsDecay>(new IsobarLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
      else if (usedSystem=="Tensor"){
	if(motherParticle->twoJ()==2 && motherParticle->theParity()==-1 && motherParticle->theCParity()==-1){
	  if (daughterParticles[0]->name()=="photon") tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorPsiToGamXDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
	  else if (daughterParticles[1]->name()=="photon") tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorPsiToGamXDecay(motherParticle, daughterParticles[1], daughterParticles[0], _channelID));
	  else tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
	}
	else tmpDec= std::shared_ptr<AbsDecay>(new IsobarTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], _channelID));
      }
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }

    }

    else if(daughterParticles.size()==3){
      if (usedSystem=="Cano") tmpDec= std::shared_ptr<AbsDecay>(new OmegaTo3PiLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], _channelID));
      else if (usedSystem=="Tensor") tmpDec= std::shared_ptr<AbsDecay>(new OmegaTo3PiTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], _channelID));
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }
    }

    else {
      Alert << "Decay\t" << (*itStr) << "\tnot supported!!!" ;
    }

    if(!useIsospin) tmpDec->disableIsospin();
    tmpDec->extractStates();
    _absDecList->addDecay(tmpDec);
    Info << "added decay " << tmpDec->name() << " to decay list" << endmsg;
  }


  //produced particle pairs
  std::vector<std::string> productionSystem = _theParser->productionSystem();

  for ( itStr = productionSystem.begin(); itStr != productionSystem.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string firstParticleStr;
    stringStr >> firstParticleStr;
    std::cout << "first particle:\t" << firstParticleStr << std::endl;

    std::string secondParticleStr;
    stringStr >> secondParticleStr;
    std::cout << "second particle:\t" << secondParticleStr << std::endl;

    Particle* firstParticle = GlobalEnv::instance()->particleTable()->particle(firstParticleStr);
    if( 0==firstParticle){
      Alert << "particle\t" << firstParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    Particle* secondParticle = GlobalEnv::instance()->particleTable()->particle(secondParticleStr);
    if( 0==secondParticle){
      Alert << "particle\t" << secondParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    _producedParticlePairs.push_back(make_pair(firstParticle, secondParticle));

  }

  //set prefactor
  std::vector<std::string> preFactorVec = _theParser->preFactor();
  for ( itStr = preFactorVec.begin(); itStr != preFactorVec.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string ampName;
    stringStr >> ampName;

    std::string valueStr;
    stringStr >> valueStr;

    double currentValue;
    currentValue=atof(valueStr.c_str());

    _preFactorMap[ampName]=currentValue;
  }


  //fill vector histMassSystems
  std::vector<std::string> theHistMassNames=_theParser->histMassNames();
  for ( itStr = theHistMassNames.begin(); itStr != theHistMassNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::vector<std::string> currentStringVec;
    while(stringStr >> tmpName){
      currentStringVec.push_back(tmpName);
    }
    _histMassSystems.push_back(currentStringVec);
  }

  //mass range cuts
  std::vector<std::string> massRangeCuts=_theParser->massRangeCuts();
  if(massRangeCuts.size()>0) _useMassRange=true;
  for ( itStr = massRangeCuts.begin(); itStr != massRangeCuts.end(); ++itStr){
    std::string currentString=*itStr;
    std::shared_ptr<MassRangeCut> currentMassRangeCut(new MassRangeCut( currentString, _finalStateParticles));
    _massRangeCuts.push_back(currentMassRangeCut);
  }  

 // hist angles
  std::vector<std::string> theHistAngleNames=_theParser->histAngleNames();
//fill vector histMassSystems
  for ( itStr = theHistAngleNames.begin(); itStr != theHistAngleNames.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::vector<std::string> currentStringDecVec;
    std::vector<std::string> currentStringDecVec2;
    std::vector<std::string> currentStringMotherVec;
    bool isDecParticle=true;
    bool isFirstDecParticle=true;
    short nBodyDecay=2;
    while(stringStr >> tmpName){
      if(tmpName=="from") {
	isDecParticle=false;
	continue;
      }
      else if(tmpName=="and") {
	isFirstDecParticle=false;
	nBodyDecay=3;
	continue;
      }
      if(isDecParticle && isFirstDecParticle) currentStringDecVec.push_back(tmpName);
      else if(isDecParticle && !isFirstDecParticle) currentStringDecVec2.push_back(tmpName);
      else currentStringMotherVec.push_back(tmpName);
    }
    std::shared_ptr<angleHistData> currentAngleHistData(new angleHistData(currentStringMotherVec, currentStringDecVec, currentStringDecVec2, nBodyDecay));
    _angleHistDataVec.push_back(currentAngleHistData);
  }

  // 2D hists
  std::vector<std::string> theHistAngleNames2D=_theParser->histAngleNames2D();
  for ( itStr = theHistAngleNames2D.begin(); itStr != theHistAngleNames2D.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::vector<std::string> currentStringDecVec;
    std::vector<std::string> currentStringMotherVec;
    std::vector<std::string> currentStringDecVec_2;
    std::vector<std::string> currentStringMotherVec_2;
    bool isDecParticle=true;
    bool isFirstAxes=true;
    while(stringStr >> tmpName){
      if(tmpName=="from") {
        isDecParticle=false;
	std::cout << tmpName << std::endl;
        continue;
      }
      else if(tmpName=="versus") {
        isFirstAxes=false;
	isDecParticle=true;
        continue;
      }
      if(isDecParticle && isFirstAxes) currentStringDecVec.push_back(tmpName);
      else if (!isDecParticle && isFirstAxes) currentStringMotherVec.push_back(tmpName);
      else if(isDecParticle && !isFirstAxes) currentStringDecVec_2.push_back(tmpName);
      else currentStringMotherVec_2.push_back(tmpName);
    }
    std::shared_ptr<angleHistData2D> currentAngleHistData2D(new angleHistData2D(currentStringMotherVec, currentStringDecVec, currentStringMotherVec_2, currentStringDecVec_2));
    _angleHistDataVec2D.push_back(currentAngleHistData2D);
  }

 // calculate contributions
  std::vector<std::string> theCalcContribution=_theParser->calcContribution();
  for ( itStr = theCalcContribution.begin(); itStr != theCalcContribution.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;
    std::string contribName;
    std::vector<std::string> currentStringZeroAmp;
    bool isContribName = true;
    while(stringStr >> tmpName){
      if(tmpName=="withZeroAmp") {
      	isContribName = false;
      	continue;
      }
      if(isContribName) contribName = tmpName;
      else currentStringZeroAmp.push_back(tmpName);
    }
    std::shared_ptr<calcContributionData> currentCalcContributionData(new calcContributionData(contribName, currentStringZeroAmp));
    _calcContributionDataVec.push_back(currentCalcContributionData);
  }


}
