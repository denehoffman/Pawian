#include <getopt.h>
#include <fstream>

#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/IsobarDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"

IsobarDecayList* IsobarDecayList::_instance=0;

IsobarDecayList* IsobarDecayList::instance()
{
  if (0==_instance) _instance = new IsobarDecayList();
  return _instance;
}

IsobarDecayList::~IsobarDecayList(){
}

void IsobarDecayList::addDecay(boost::shared_ptr<IsobarDecay> theIsoDec){
  const std::string key = theIsoDec->motherPart()->name();

  //  std::map<std::string, boost::shared_ptr<IsobarDecay> >::iterator it;
  //  it=_isoDecList.find(key);
  //  if (it != _isoDecList.end()){
    _isoDecList[theIsoDec->motherPart()->name()]=theIsoDec;
    //  }
    //  else{
    //    Warning << "decay of the particle " << key << " already exists\n"
    //	    << "new decay cannot be added!!!" << endmsg;  
    //  }
}

boost::shared_ptr<IsobarDecay> IsobarDecayList::decay(Particle* mother){
   
  boost::shared_ptr<IsobarDecay> result;

  std::string key = mother->name();

  //   std::map<std::string, boost::shared_ptr<IsobarDecay> >::iterator it;
  //   it=_isoDecList.find(key);
  //   if (it!=_isoDecList.end()){
    result=_isoDecList[key];
    //     }
    //     else{
    // Warning << "decay of the particle " << key << " does not exist"
    //	    << endmsg;  
    // }

  return result; 
}
