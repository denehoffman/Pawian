//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include "PwaUtils/PhpGenDynamics.hh"
#include "PwaDynamics/BreitWignerFunction.hh"
#include "Particle/Particle.hh"

PhpGenDynamics::PhpGenDynamics(std::string& description, std::vector<Particle*>& finalStateParticles) :
  _dynName("")
  ,_mass0(0.)
  ,_width0(0.)
  ,_massDaughter1(0.)
  ,_massDaughter2(0.)
  ,_maxWeight(0.) 
{
  // required order for the string: massMin massMax particleName1 particleName2 ...
  std::vector<Particle*>::iterator itParticle;
  int counter=0;
  std::stringstream stringStrDescription;
  stringStrDescription << description;
  std::string tmpDescription;
  while(stringStrDescription >> tmpDescription){
    InfoMsg <<"\nPhpGenDynamics description content:\t" << tmpDescription << endmsg;
    if (counter==0) _dynName= tmpDescription;
    else if(counter==1) _mass0=atof(tmpDescription.c_str());
    else if(counter==2) _width0=atof(tmpDescription.c_str());
    else{
      unsigned int index=0;
      for (itParticle=finalStateParticles.begin(); itParticle != finalStateParticles.end(); ++itParticle){
	//	InfoMsg << "tmpNameMassRange: " << tmpNameMassRange << "\t(*itParticle)->name(): " << (*itParticle)->name() <<endmsg;
	if( tmpDescription==(*itParticle)->name()){
	  //	  InfoMsg << "particle " << tmpNameMassRange  << " added to list of mass range cuts " <<endmsg;
	  _decParticles.push_back(*itParticle);
	  _iDparticles.push_back(index);
	  break;	  
	}
	++index;
      } 
    }
    ++counter;
  }

  if(_mass0<0. || _width0 <0.){
    Alert << "_mass0<0. or _width0 < 0 !!!!!" << endmsg;
    exit(0);
  }

  if(_decParticles.size()<2){
    Alert << "Al least 2 decay particles are required!" 
	  << "\n here _decParticles.size()= " <<  _decParticles.size() << endmsg;
    exit(0);
  }
  
  if (_dynName=="BreitWigner"){
    _maxWeight=eval(_mass0);
  }
  else{
    Alert << "dynamics: " << _dynName << " is not supported!!!" << endmsg;
    exit(0);
  } 


  InfoMsg << "mass cut on inv mass of the final state particles:" <<endmsg;
  for (itParticle=_decParticles.begin(); itParticle!=_decParticles.end(); ++itParticle){
    InfoMsg << (*itParticle)->name() << endmsg;
  }

  InfoMsg << "\tdynamics in php generation:\t" << _dynName 
	  << endmsg;
}

PhpGenDynamics::~PhpGenDynamics()
{
}

bool PhpGenDynamics::isDecayParticle(Particle* compParticle){
  bool result=false;
  std::vector<Particle*>::iterator it;
  for(it=_decParticles.begin(); it!=_decParticles.end(); ++it){
    if ((*it)->sameName( *compParticle )){
      result=true;
      break;
    }   
  } 
  return result;
}

double PhpGenDynamics::eval(double currentMass){
  return norm(BreitWignerFunction::NonRel(currentMass,_mass0, _width0));
} 




