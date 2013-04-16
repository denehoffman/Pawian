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

// AbsEnv class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <stdlib.h>

#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "PwaUtils/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/PdtParser.hh"


AbsEnv::AbsEnv() :
  _alreadySetUp(false)
  , _noFinalStateParticles(0)
  ,_absDecList(new AbsDecayList())
  ,_prodDecList(new AbsDecayList())
  ,_useMassRange(false)
  ,_massMin(0.)
  ,_massMax(100.)
  ,_theParser(0)
{
}

AbsEnv::~AbsEnv(){
}

void AbsEnv::setup(ParserBase* theParser){
  if(_alreadySetUp){
    Alert << " AbsEnv already set up" << endmsg;
    exit(1);
  }
  
  _alreadySetUp=true;

  _theParser=theParser;
  // common options
  _outputFileNameSuffix = theParser->outputFileNameSuffix();
  _serializationFileName = theParser->serializationFile();

  // pdtTable
  PdtParser pdtParser;
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  std::string pdtFileRelPath=theParser->pdgTableFile(); 
  std::string pdtFile(theSourcePath+pdtFileRelPath);
  _particleTable = new ParticleTable;
  
  if (!pdtParser.parse(pdtFile, *_particleTable)) {
    Alert << "can not parse particle table " << pdtFile << endmsg;
    exit(1);
  }

  // cloned particles
  const std::vector<std::string> cloneParticle=theParser->cloneParticle();
  std::vector<std::string>::const_iterator itcP;

  for ( itcP = cloneParticle.begin(); itcP != cloneParticle.end(); ++itcP){
     std::istringstream particles(*itcP);
     std::string particleOld;
     std::string particleNew;
     particles >> particleOld >> particleNew;

     _particleTable->clone(particleNew, particleOld);
  }


  //final state particles
  const std::vector<std::string> finalStateParticleStr=theParser->finalStateParticles();
  
  std::vector<std::string>::const_iterator itStr;
  for ( itStr = finalStateParticleStr.begin(); itStr != finalStateParticleStr.end(); ++itStr){
    Particle* currentParticle = _particleTable->particle(*itStr);
    _finalStateParticles.push_back(currentParticle);
  }

  _noFinalStateParticles= (int) _finalStateParticles.size();

  //decays

  std::vector<std::string> decaySystem= theParser->decaySystem();
  for ( itStr = decaySystem.begin(); itStr != decaySystem.end(); ++itStr){

    Particle* motherParticle =0;
    std::vector<Particle*> daughterParticles;

    std::stringstream stringStr;
    stringStr << (*itStr);

    std::string tmpName;

    bool isDecParticle=false;
    bool firstArgument=true;
    std::string usedSystem("");

    while(stringStr >> tmpName){
      if(firstArgument){
	usedSystem=tmpName;
	firstArgument=false;
	continue;
      }
      if(tmpName=="To") {
        isDecParticle=true;
        continue;
      }
      if(isDecParticle){
  	daughterParticles.push_back(_particleTable->particle(tmpName));
      }
      else{
  	motherParticle = _particleTable->particle(tmpName);
      }
    }
    boost::shared_ptr<AbsDecay> tmpDec;
    if(daughterParticles.size()==2){
      if (usedSystem=="Heli") tmpDec= boost::shared_ptr<AbsDecay>(new IsobarHeliDecay(motherParticle, daughterParticles[0], daughterParticles[1], this));
      else if (usedSystem=="Cano")  tmpDec= boost::shared_ptr<AbsDecay>(new IsobarLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], this));
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }

    }

    else if(daughterParticles.size()==3){
      if (usedSystem=="Cano") tmpDec= boost::shared_ptr<AbsDecay>(new OmegaTo3PiLSDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], this));
      else if (usedSystem=="Tensor") tmpDec= boost::shared_ptr<AbsDecay>(new OmegaTo3PiTensorDecay(motherParticle, daughterParticles[0], daughterParticles[1], daughterParticles[2], this));
      else {
	Alert << "used decay system\t" << usedSystem << "\tnot supported!!!\n" << endmsg;
	exit(1);
      }
    }
 
    else {
      Alert << "Decay\t" << (*itStr) << "\tnot supported!!!" ; 
    }

     _absDecList->addDecay(tmpDec);
  }


  //produced particle pairs
  std::vector<std::string> productionSystem = theParser->productionSystem();

  for ( itStr = productionSystem.begin(); itStr != productionSystem.end(); ++itStr){
    std::stringstream stringStr;
    stringStr << (*itStr);
    std::string firstParticleStr;
    stringStr >> firstParticleStr;
    std::cout << "first particle:\t" << firstParticleStr << std::endl;

    std::string secondParticleStr;
    stringStr >> secondParticleStr;
    std::cout << "second particle:\t" << secondParticleStr << std::endl;

    Particle* firstParticle = _particleTable->particle(firstParticleStr);
    if( 0==firstParticle){
      Alert << "particle\t" << firstParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    Particle* secondParticle = _particleTable->particle(secondParticleStr);
    if( 0==secondParticle){
      Alert << "particle\t" << secondParticleStr << "\tdoes not exist in pdtTable" << endmsg;
      exit(1);
    }

    _producedParticlePairs.push_back(make_pair(firstParticle, secondParticle));

  }


  //fill vector histMassSystems
  std::vector<std::string> theHistMassNames=theParser->histMassNames();
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

  //mass range
  int counter=0;
  std::string massRangeStr=theParser->massRange();
  if(massRangeStr.size()>0) _useMassRange=true;

  std::stringstream stringStrMassRange;
  stringStrMassRange << massRangeStr;
  std::string tmpNameMassRange;

  while(stringStrMassRange >> tmpNameMassRange){
    Info <<"\nMassRangeCont:\t" << tmpNameMassRange << endmsg;
    if (counter==0) _massMin=atof(tmpNameMassRange.c_str());
    else if(counter==1) _massMax=atof(tmpNameMassRange.c_str());
    else{
      //find index
      for(size_t idex=0; idex<_finalStateParticles.size(); ++idex){
	Particle* currentParticle=_finalStateParticles[idex];
	if(currentParticle->name() == tmpNameMassRange){
	  _particleIndicesMassRange.push_back(idex);
	  Info << "\nFound particle\t" << currentParticle->name() << "\t index:\t" << idex << endmsg; 
	  break;
	}
      }
    }
    counter++;
  }
}

