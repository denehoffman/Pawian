#include <getopt.h>
#include <fstream>

#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


IsobarDecay::IsobarDecay(Particle* mother, Particle* daughter1, Particle* daughter2) :
  _mother(mother)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_motherJPCPtr(getJPCPtr(mother)) 

  ,_daughter1JPCPtr(getJPCPtr(daughter1))
  ,_daughter2JPCPtr(getJPCPtr(daughter2))
  ,_name(mother->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
{
  validJPCLS( _motherJPCPtr, _daughter1JPCPtr, _daughter2JPCPtr, _JPCLSDecAmps);
  _isoDecDaughter1=pbarpEnv::instance()->decayList()->decay(_daughter1);
  _isoDecDaughter2=pbarpEnv::instance()->decayList()->decay(_daughter2);
}

IsobarDecay::IsobarDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, std::string motherName) :
  _mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_motherJPCPtr(motherJPCPtr)
  ,_daughter1JPCPtr(getJPCPtr(daughter1))
  ,_daughter2JPCPtr(getJPCPtr(daughter2))
  ,_name(motherName+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
{
  validJPCLS( _motherJPCPtr, _daughter1JPCPtr, _daughter2JPCPtr, _JPCLSDecAmps);
  _isoDecDaughter1=pbarpEnv::instance()->decayList()->decay(_daughter1);
  _isoDecDaughter2=pbarpEnv::instance()->decayList()->decay(_daughter2);
}

IsobarDecay::~IsobarDecay(){
}

void IsobarDecay::print(std::ostream& os) const{
  os << "\nJPCLS amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator it;
  for (it = _JPCLSDecAmps.begin(); it!= _JPCLSDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }
  
  if(0!=_isoDecDaughter1){
    os << "with further decay:";
    _isoDecDaughter1->print(os);
  }
  if(0!=_isoDecDaughter2){
    os << "with further decay:";
    _isoDecDaughter2->print(os);
  }  

  os << "\n";
}
