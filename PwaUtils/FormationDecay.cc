//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// FormationDecay class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/FormationDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

FormationDecay::FormationDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, ChannelID channelID, std::string motherName, std::string typeName) :
  AbsDecay(motherIGJPCPtr, daughter1, motherName, channelID, typeName)
{
}

FormationDecay::~FormationDecay(){
}

void FormationDecay::print(std::ostream& os) const{
  AbsDecay::print(os);
  os << "\n";
}

void FormationDecay::extractStates(){
}

void FormationDecay::setDecayLevel(decLevel theLevel){
  if (_decLevel==decayLevel::noLevel){
    _decLevel=theLevel;
    if(_decLevel==decayLevel::isProdAmp) _hasMotherPart=false;
  }
  else if(_decLevel != theLevel) _decLevel=decayLevel::severalLevels;
  else return;
  InfoMsg << name() << " set decay level to " << _decLevel << endmsg; 
}

void FormationDecay::setDecayLevelTree(decLevel theLevel, std::shared_ptr<AbsDecay> motherDecPtr, std::shared_ptr<AbsDecay> prodDecPtr){
  setDecayLevel(theLevel);
  if(prodDecPtr->whichDecayLevel() != AbsDecay::decayLevel::isProdAmp){
    Alert << "prodDecPtr with the name " << prodDecPtr->name() << " is not a production amplitude!!!" << endmsg;
    exit(1);
  }

  //check if prod decay already exists
  const std::string key = prodDecPtr->name();
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for(it=_prodAmpRefList.begin(); it!=_prodAmpRefList.end();++it){
    if (key==(*it)->name()){
      //prod decay already exists
      Alert << "production amplitude " << key << " aleady exists in the reference list!!!" << endmsg;
      Alert << "It is not allowed to refer to the same production amplitude more than once!!!!" << endmsg;
      exit(1);
    }
  }

  _prodAmpRefList.push_back(prodDecPtr);

  const std::string keyMotherDec = motherDecPtr->name();
  for(it=_motherAmpRefList.begin(); it!=_motherAmpRefList.end();++it){
    if (key==(*it)->name()){
      //prod decay already exists
      Alert << "mother amplitude " << keyMotherDec << " aleady exists in the reference list!!!" << endmsg;
      Alert << "It is not allowed to refer to the same mother amplitude more than once!!!!" << endmsg;
      exit(1);
    }
  }

  if(whichDecayLevel() != AbsDecay::decayLevel::isProdAmp) _motherAmpRefList.push_back(motherDecPtr);

  _absDecDaughter1->setDecayLevelTree(decayLevel(theLevel+1), shared_from_this(), prodDecPtr);
}

void FormationDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData){

  std::vector<Particle*>::iterator itP;
  std::map<std::string, Vector4<double> >::iterator itMap;
  
  Vector4<double> mother4Vec(0.,0.,0.,0.);
  
  //fill mother4Vec
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    mother4Vec+=itMap->second;
  }

  _absDecDaughter1->fillWignerDs(fsMap, mother4Vec, evtData); 
}
