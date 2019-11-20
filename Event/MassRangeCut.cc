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

#include "Event/MassRangeCut.hh"
#include "Particle/Particle.hh"

MassRangeCut::MassRangeCut(std::string& rangeAndParticleNames, 
			   std::vector<Particle*>& finalStateParticles, 
			   bool anticut) : 
  _massMin(0.),
  _massMax(0.),
  _isAntiCut(anticut)
{
  // required order for the string: massMin massMax particleName1 particleName2 ...
  std::vector<Particle*>::iterator itParticle;
  int counter=0;
  std::stringstream stringStrMassRange;
  stringStrMassRange << rangeAndParticleNames;
  std::string tmpNameMassRange;
  while(stringStrMassRange >> tmpNameMassRange){
    InfoMsg <<"\nMassRangeCont:\t" << tmpNameMassRange << endmsg;
    if (counter==0) _massMin=atof(tmpNameMassRange.c_str());
    else if(counter==1) _massMax=atof(tmpNameMassRange.c_str());
    else{
      unsigned int index=0;
      for (itParticle=finalStateParticles.begin(); itParticle != finalStateParticles.end(); ++itParticle) {
	//	InfoMsg << "tmpNameMassRange: " << tmpNameMassRange << "\t(*itParticle)->name(): " << (*itParticle)->name() <<endmsg;
	if( tmpNameMassRange==(*itParticle)->name()){
	  //	  InfoMsg << "particle " << tmpNameMassRange  << " added to list of mass range cuts " <<endmsg;
	  _particlesMassRange.push_back(*itParticle);
	  _iDparticlesMassRange.push_back(index);
	  break;	  
	}
	++index;
      } 
    }
    ++counter;
  }
  if(_massMin<0. || _massMin>_massMax){
    Alert << "massMin<0. or min>massMax !!!!!" << endmsg;
    exit(0);
  }
  if(_particlesMassRange.size()<2){
    Alert << "It is required to cut on an invariant mass of at least 2 final state particles!" 
	  << "\n here particlesMassRange.size()= " <<  _particlesMassRange.size() << endmsg;
    exit(0);
  }

  InfoMsg << "mass cut on inv mass of the final state particles:" <<endmsg;
  for (itParticle=_particlesMassRange.begin(); itParticle!=_particlesMassRange.end(); ++itParticle){
    InfoMsg << (*itParticle)->name() << endmsg;
  }

  InfoMsg << "\tmassMin = " << _massMin 
	  << "\tmassmax = " << _massMax
	  << endmsg;
}

MassRangeCut::~MassRangeCut()
{
}

bool MassRangeCut::isMassRangeParticle(Particle* compParticle){
  bool result=false;
  std::vector<Particle*>::iterator it;
  for(it=_particlesMassRange.begin(); it!=_particlesMassRange.end(); ++it){
    if ((*it)->sameName( *compParticle )){
      result=true;
      break;
    }   
  } 
  return result;
} 




