#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "Examples/pbarp/IsobarDecay.hh"
#include "Examples/pbarp/IsobarDecayList.hh"
#include "Examples/pbarp/pbarpEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"

IsobarDecay::IsobarDecay(Particle* mother, Particle* daughter1, Particle* daughter2) :
  _mother(mother)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(true)
  ,_motherJPCPtr(getJPCPtr(mother)) 
  ,_daughter1JPCPtr(getJPCPtr(daughter1))
  ,_daughter2JPCPtr(getJPCPtr(daughter2))
  ,_name(mother->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
{
  validJPCLS( _motherJPCPtr, _daughter1JPCPtr, _daughter2JPCPtr, _JPCLSDecAmps);
  _isoDecDaughter1=pbarpEnv::instance()->decayList()->decay(_daughter1);
  if(0 != _isoDecDaughter1){
    _daughter1IsStable=false;
    std::vector<Particle*> fsParticlesDaughter1=_isoDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), fsParticlesDaughter1.begin(), fsParticlesDaughter1.end());
  }
  else _finalStateParticles.push_back(daughter1);
  _isoDecDaughter2=pbarpEnv::instance()->decayList()->decay(_daughter2);

  if(0 != _isoDecDaughter2){
    _daughter2IsStable=false;
    _finalStateParticlesDaughter2=_isoDecDaughter2->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter2.begin(), _finalStateParticlesDaughter2.end());
  }
  else{
    _finalStateParticles.push_back(daughter2);
    _finalStateParticlesDaughter2.push_back(daughter2);
  }

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);
  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+_motherJPCPtr->name()+FunctionUtils::particleListName(_finalStateParticles);

}

IsobarDecay::IsobarDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, std::string motherName) :
  _mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(false)
  ,_motherJPCPtr(motherJPCPtr)
  ,_daughter1JPCPtr(getJPCPtr(daughter1))
  ,_daughter2JPCPtr(getJPCPtr(daughter2))
  ,_name(motherName+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
{
  validJPCLS( _motherJPCPtr, _daughter1JPCPtr, _daughter2JPCPtr, _JPCLSDecAmps);
  _isoDecDaughter1=pbarpEnv::instance()->decayList()->decay(_daughter1);


  if(0 != _isoDecDaughter1){
    _daughter1IsStable=false;
    std::vector<Particle*> fsParticlesDaughter1=_isoDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), fsParticlesDaughter1.begin(), fsParticlesDaughter1.end());
  }
  else _finalStateParticles.push_back(daughter1);
  _isoDecDaughter2=pbarpEnv::instance()->decayList()->decay(_daughter2);

  if(0 != _isoDecDaughter2){
    _daughter2IsStable=false;
    _finalStateParticlesDaughter2=_isoDecDaughter2->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter2.begin(), _finalStateParticlesDaughter2.end());
  }
  else{
    _finalStateParticles.push_back(daughter2);
    _finalStateParticlesDaughter2.push_back(daughter2);
  }
  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+motherName;
}

IsobarDecay::~IsobarDecay(){
}

void IsobarDecay::fillWignerDs(std::map<std::string, Vector4<float> >& fsMap, EvtDataNew* evtData){
  if (!_daughter1IsStable) _isoDecDaughter1->fillWignerDs(fsMap, evtData);
  if (!_daughter2IsStable) _isoDecDaughter2->fillWignerDs(fsMap, evtData);
  
  Vector4<float> all4Vec(0.,0.,0.,0.);
  Vector4<float> mother4Vec(0.,0.,0.,0.);
  Vector4<float> daughter2_4Vec(0.,0.,0.,0.);

  //fill all4Vec
  std::map<std::string, Vector4<float>, pawian::Collection::PtrLess>::iterator itMap;
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill mother4Vec
  std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    mother4Vec+=itMap->second;
  }

  //fill mother4Vec
  for(itP=_finalStateParticlesDaughter2.begin(); itP!=_finalStateParticlesDaughter2.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    daughter2_4Vec+=itMap->second;
  }

  Vector4<float> daughter2HelMother(0.,0.,0.,0.);
  if(_hasMotherPart){
    daughter2HelMother=helicityVec(all4Vec, mother4Vec, daughter2_4Vec);
  }
  else{
    daughter2HelMother=daughter2_4Vec;
    daughter2HelMother.Boost(mother4Vec);
  }

  Spin spinMother=_motherJPCPtr->J;
  Spin spinDaughter1=_daughter1JPCPtr->J;
  Spin spinDaughter2=_daughter2JPCPtr->J;
  Spin lam12Max=spinDaughter1+spinDaughter2;
  if(lam12Max>spinMother) lam12Max=spinMother;

  Spin lamMotherMax=spinMother;
  if (!_hasMotherPart && spinMother>1) lamMotherMax=1;

 
  for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
    for (Spin lam12=-lam12Max; lam12<=lam12Max; ++lam12){
      double thePhi=0.;
      if(_hasMotherPart) thePhi=daughter2HelMother.Phi();
      evtData->WignerDsString[_wignerDKey][spinMother][lamMother][lam12]=Wigner_D(thePhi,daughter2HelMother.Theta(),0,spinMother,lamMother,lam12);
    }
  }
}

void IsobarDecay::print(std::ostream& os) const{
  os << "\nJPCLS amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator it;
  for (it = _JPCLSDecAmps.begin(); it!= _JPCLSDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }
  
  if(!_daughter1IsStable){
    os << "with further decay:";
    _isoDecDaughter1->print(os);
  }
  if(!_daughter2IsStable){
    os << "with further decay:";
    _isoDecDaughter2->print(os);
  }  

  os << "\n";
}
