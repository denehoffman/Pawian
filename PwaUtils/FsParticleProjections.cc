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

// FsParticleProjections class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/FsParticleProjections.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

FsParticleProjections::FsParticleProjections(std::vector<Particle*>& fsParticles) :
  _fsParticles(fsParticles)
{
  _fsParticleNames.resize(_fsParticles.size());
  _fsSpins.resize(_fsParticles.size());

  unsigned int noOfProjections=1;
  std::vector< std::vector<Spin> > fsSpinProjectionsInverse;
  fsSpinProjectionsInverse.resize(_fsParticles.size());
  std::vector<unsigned int > spinMaxArrays;
  spinMaxArrays.resize(_fsParticles.size());
  std::vector<unsigned int > currentSpinProjections;
  currentSpinProjections.resize(_fsParticles.size());

  for(unsigned int i=0; i<_fsParticles.size(); ++i){
    _fsParticleNames[i]=_fsParticles.at(i)->name();
    double theSpin=_fsParticles.at(i)->twoJ()/2.;
    Spin currentSpin(theSpin);
    _fsSpins[i]=currentSpin;
    std::vector<Spin> currentProjections;
    Spin currentProjection(-currentSpin);
    while (currentProjection <= currentSpin){
      //      if(currentProjection==Spin(0) && _fsParticles.at(i)->mass() < 1.e-10) continue;
      if(currentProjection != Spin(0)) currentProjections.push_back(currentProjection);
      else if (_fsParticles.at(i)->mass() > 1.e-9) currentProjections.push_back(currentProjection);
      currentProjection++;
    }
    fsSpinProjectionsInverse[i]=currentProjections;
    noOfProjections*=currentProjections.size();
    spinMaxArrays[i]=currentProjections.size();
    
    currentSpinProjections[i]=0;    
  }

  Info << "\nThe spins of the final state particles are:" << endmsg;
  for(unsigned int i=0; i<_fsParticles.size(); ++i){
    Info << _fsParticleNames.at(i) << " : " << _fsSpins.at(i) << endmsg; 
  }

  unsigned int currentPIndex=0;
  for(unsigned int i=0; i<noOfProjections; ++i){

    std::vector<Spin> currentSpinHelicities;
    currentSpinHelicities.resize(_fsParticles.size());
    Info <<"\ncurrent projection is:" << endmsg;
    for (unsigned int ii=0; ii<_fsParticles.size(); ++ii){
      Info << currentSpinProjections.at(ii) << "\tlambda: " << fsSpinProjectionsInverse.at(ii).at(currentSpinProjections.at(ii)) << endmsg;
      currentSpinHelicities[ii]=fsSpinProjectionsInverse.at(ii).at(currentSpinProjections.at(ii));
    }
    _fsSpinProjections.push_back(currentSpinHelicities);
    if(i<noOfProjections-1) increaseCurrentSpinIndex(currentSpinProjections, spinMaxArrays, currentPIndex);
  }
}


FsParticleProjections::~FsParticleProjections()
{
}

bool FsParticleProjections::increaseCurrentSpinIndex(std::vector<unsigned int >& currentSpinProjections, std::vector<unsigned int >& spinMaxArrays, unsigned int& currentPIndex){
  bool result = false;
  for(unsigned int id=0; id<=currentPIndex; ++id){
    if ( currentSpinProjections.at(id) < (spinMaxArrays.at(currentPIndex)-1) ){
      currentSpinProjections.at(id)++;
      if ( currentSpinProjections.at(id) < (spinMaxArrays.at(currentPIndex)-1) ) currentPIndex=id;
      else{
	for(unsigned int idReset=0; idReset<id; ++idReset){
	  currentSpinProjections.at(idReset) = 0;
	  currentPIndex=0;
	}
      }
      return true;
    }
  }

  currentPIndex++;
  while (currentSpinProjections.at(currentPIndex) >= (spinMaxArrays.at(currentPIndex)-1))  currentPIndex++;
  currentSpinProjections.at(currentPIndex)++;
  for(unsigned int idReset=0; idReset< currentPIndex; ++idReset){
    currentSpinProjections.at(idReset) = 0;
  } 
  Info << "currentPIndex: " << currentPIndex << endmsg;
  result = true; 
  return result; 
}

unsigned int FsParticleProjections::fsParticleId(std::string& name){
  unsigned int result=_fsParticleNames.size()+1;
  for(unsigned int id=0; id<_fsParticleNames.size(); ++id){
    Info << "_fsParticleNames.at(" << id <<"): " << _fsParticleNames.at(id) << endmsg;
    if(name == _fsParticleNames.at(id)){
      result=id;
      break;
    } 
  }
  if(result>_fsParticleNames.size()){
    Alert << "particle with name " << name << " does not exist!!!" << endmsg;
    exit(1);
  }
  return result;
}



