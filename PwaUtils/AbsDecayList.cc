// AbsDecayList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


AbsDecayList::AbsDecayList(){
}

AbsDecayList::~AbsDecayList(){
}

void AbsDecayList::addDecay(boost::shared_ptr<AbsDecay> theIsoDec){

//check if this decay already exists
  bool alreadyThere=false;  
  const std::string key = theIsoDec->name();
  std::vector<boost::shared_ptr<AbsDecay> >::iterator it;
  for(it=_absDecList.begin(); it!=_absDecList.end();++it){
    if (key==(*it)->name()){
      alreadyThere=true;
      break;
    }
  }
  
  if(!alreadyThere) _absDecList.push_back(theIsoDec);
  else{
    Warning << "decay of the particle " << key << " already exists\n"
	    << "new decay cannot be added!!!" << endmsg;  
  }
  
}

boost::shared_ptr<AbsDecay> AbsDecayList::decay(Particle* mother){
   
  boost::shared_ptr<AbsDecay> result;

  const std::string key = mother->name();

  std::vector<boost::shared_ptr<AbsDecay> >::iterator it;
  for (it=_absDecList.begin(); it!=_absDecList.end(); ++it){
    if (key==(*it)->motherPart()->name()){
      result=(*it);
      break;
    }

  }
  return result; 
}

boost::shared_ptr<AbsDecay> AbsDecayList::decay(const std::string& name){
  boost::shared_ptr<AbsDecay> result;

  std::vector<boost::shared_ptr<AbsDecay> >::iterator it;
  for (it=_absDecList.begin(); it!=_absDecList.end(); ++it){
    if (name==(*it)->name()){
      result=(*it);
      break;
    }
    
  }

  return result; 

}

void AbsDecayList::replaceSuffix(const std::string& oldPart, const std::string& newPart){

  std::vector<boost::shared_ptr<AbsDecay> >::iterator it;
  for (it= _absDecList.begin(); it!=_absDecList.end(); ++it){
    std::string theSuffix= (*it)->fitParSuffix();
    std::cout << "theSuffix:\t" << theSuffix << std::endl;
    boost::replace_all(theSuffix, oldPart, newPart);
    (*it)->setFitParSuffix(theSuffix); 
  }
}

void AbsDecayList::replaceMassKey(const std::string& oldPart, const std::string& newPart){
  std::vector<boost::shared_ptr<AbsDecay> >::iterator it;
  for (it= _absDecList.begin(); it!=_absDecList.end(); ++it){
    if(oldPart== (*it)->massParKey()){
      (*it)->setMassParKey(newPart);
    } 
  }
}
