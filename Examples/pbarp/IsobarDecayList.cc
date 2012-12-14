#include <getopt.h>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


IsobarDecayList::IsobarDecayList(){
}

IsobarDecayList::~IsobarDecayList(){
}

void IsobarDecayList::addDecay(boost::shared_ptr<IsobarDecay> theIsoDec){
  //  const std::string key = theIsoDec->motherPart()->name();

  // std::map<const std::string, boost::shared_ptr<IsobarDecay> >::iterator it;
  // it=_isoDecList.find(key);
  // if (it != _isoDecList.end()){
  _isoDecList.push_back(theIsoDec);
  // }
  // else{
  //   Warning << "decay of the particle " << key << " already exists\n"
  //   	    << "new decay cannot be added!!!" << endmsg;  
  // }
}

boost::shared_ptr<IsobarDecay> IsobarDecayList::decay(Particle* mother){
   
  boost::shared_ptr<IsobarDecay> result;

  const std::string key = mother->name();

  std::vector<boost::shared_ptr<IsobarDecay> >::iterator it;
  for (it=_isoDecList.begin(); it!=_isoDecList.end(); ++it){
    if (key==(*it)->motherPart()->name()){
      result=(*it);
      break;
    }

  }
  return result; 
}

boost::shared_ptr<IsobarDecay> IsobarDecayList::decay(const std::string& name){
  boost::shared_ptr<IsobarDecay> result;

  std::vector<boost::shared_ptr<IsobarDecay> >::iterator it;
  for (it=_isoDecList.begin(); it!=_isoDecList.end(); ++it){
    if (name==(*it)->name()){
      result=(*it);
      break;
    }
    
  }

  return result; 

}

void IsobarDecayList::replaceSuffix(const std::string& oldPart, const std::string& newPart){

  std::vector<boost::shared_ptr<IsobarDecay> >::iterator it;
  for (it= _isoDecList.begin(); it!=_isoDecList.end(); ++it){
    std::string theSuffix= (*it)->fitParSuffix();
    boost::replace_all(theSuffix, oldPart, newPart);
    (*it)->setFitParSuffix(theSuffix); 
  }
}
