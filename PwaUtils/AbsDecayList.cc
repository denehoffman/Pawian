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

// AbsDecayList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


AbsDecayList::AbsDecayList(){
}

AbsDecayList::~AbsDecayList(){
}

void AbsDecayList::addDecay(std::shared_ptr<AbsDecay> theIsoDec){

//check if this decay already exists
  bool alreadyThere=false;  
  const std::string key = theIsoDec->name();
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for(it=_absDecList.begin(); it!=_absDecList.end();++it){
    if (key==(*it)->name()){
      alreadyThere=true;
      break;
    }
  }
  
  if(!alreadyThere) _absDecList.push_back(theIsoDec);
  else{
    WarningMsg << "decay of the particle " << key << " already exists\n"
	    << "new decay cannot be added!!!" << endmsg;  
  }
  
}

std::shared_ptr<AbsDecay> AbsDecayList::decay(Particle* mother){
   
  std::shared_ptr<AbsDecay> result;

  const std::string key = mother->name();

  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for (it=_absDecList.begin(); it!=_absDecList.end(); ++it){
    if (key==(*it)->motherPart()->name()){
      result=(*it);
      break;
    }

  }
  return result; 
}

std::shared_ptr<AbsDecay> AbsDecayList::decay(const std::string& name){
  std::shared_ptr<AbsDecay> result;

  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for (it=_absDecList.begin(); it!=_absDecList.end(); ++it){
    if (name==(*it)->name()){
      result=(*it);
      break;
    }
    
  }

  return result; 

}

void AbsDecayList::replaceSuffix(const std::string& oldPart, const std::string& newPart){

  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for (it= _absDecList.begin(); it!=_absDecList.end(); ++it){
    std::string theSuffix= (*it)->fitParSuffix();
    std::string oldSuffix=theSuffix;
    boost::replace_all(theSuffix, oldPart, newPart);
    (*it)->setFitParSuffix(theSuffix);
    if(theSuffix != oldSuffix){
      GlobalEnv::instance()->addToStringStringMap(oldPart, newPart, GlobalEnv::instance()->alreadyReplacedSuffixMap());
      GlobalEnv::instance()->addToStringStringMap(oldSuffix, theSuffix, GlobalEnv::instance()->fitParamReplacementMap());      
    }
  }

}

void AbsDecayList::replaceMassKey(const std::string& oldPart, const std::string& newPart){
  bool noReplacement=true;
  bool alreadyReplaced=false;
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for (it= _absDecList.begin(); it!=_absDecList.end(); ++it){
    if(newPart == (*it)->massParKey()) alreadyReplaced=true;
    if(oldPart == (*it)->massParKey()){
      (*it)->setMassParKey(newPart);
      noReplacement=false;
    } 
  }

  if(noReplacement && !alreadyReplaced){
    Alert << "mass key: " << oldPart << " does not exist and cannot be replaced by the new key: " << newPart << endmsg;
    exit(1); 
  }
}

void AbsDecayList::replaceProdKey(const std::string& oldPart, const std::string& newPart){
  bool noReplacement=true;
  bool alreadyReplaced=false;
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for (it= _absDecList.begin(); it!=_absDecList.end(); ++it){
    if(newPart == (*it)->prodParKey()) alreadyReplaced=true;
    if(oldPart== (*it)->prodParKey()){
      (*it)->setProdParKey(newPart);
      noReplacement=false;
    } 
  }

  if(noReplacement && !alreadyReplaced){
    Alert << "production key: " << oldPart << " does not exist and cannot be replaced by the new key: " << newPart << endmsg;
    exit(1); 
  }
}

void AbsDecayList::print(){
  InfoMsg <<"***** the decay channels are ******** " << endmsg;
  std::vector<std::shared_ptr<AbsDecay> >::iterator it;
  for (it= _absDecList.begin(); it!=_absDecList.end(); ++it){
    InfoMsg << (*it)->name() << endmsg;
  }
}
