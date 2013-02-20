// AbsDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

AbsDecay::AbsDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv) :
  _mother(mother)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(true)
  ,_withDyn(false)
  ,_motherJPCPtr(getJPCPtr(mother)) 
  ,_daughter1JPCPtr(getJPCPtr(daughter1))
  ,_daughter2JPCPtr(getJPCPtr(daughter2))
  ,_name(mother->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
  ,_massParamKey("WoDynamics")
  ,_env(theEnv)
{
  _absDecDaughter1=_env->absDecayList()->decay(_daughter1);
  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    _finalStateParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter1.begin(), _finalStateParticlesDaughter1.end());
  }
  else{
    _finalStateParticles.push_back(daughter1);
    _finalStateParticlesDaughter1.push_back(daughter1);
  }

  _absDecDaughter2=_env->absDecayList()->decay(_daughter2);

  if(0 != _absDecDaughter2){
    _daughter2IsStable=false;
    _finalStateParticlesDaughter2=_absDecDaughter2->finalStateParticles();
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

AbsDecay::AbsDecay(boost::shared_ptr<const jpcRes> motherJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName) :
  _mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(false)
  ,_withDyn(false)
  ,_motherJPCPtr(motherJPCPtr)
  ,_daughter1JPCPtr(getJPCPtr(daughter1))
  ,_daughter2JPCPtr(getJPCPtr(daughter2))
  ,_name(motherName+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
  ,_massParamKey("WoDynamics")
  ,_env(theEnv)
{
  _absDecDaughter1=_env->absDecayList()->decay(_daughter1);


  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    std::vector<Particle*> fsParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), fsParticlesDaughter1.begin(), fsParticlesDaughter1.end());
  }
  else _finalStateParticles.push_back(daughter1);
  _absDecDaughter2=_env->absDecayList()->decay(_daughter2);

  if(0 != _absDecDaughter2){
    _daughter2IsStable=false;
    _finalStateParticlesDaughter2=_absDecDaughter2->finalStateParticles();
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

AbsDecay::~AbsDecay(){
}


void AbsDecay::enableDynamics(std::string& dynString) {
  _withDyn=true;
  if(0!=_mother) _massParamKey=_mother->name(); 
  else _massParamKey=_motherJPCPtr->name();
}

void AbsDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, EvtData* evtData){
  if (!_daughter1IsStable) _absDecDaughter1->fillWignerDs(fsMap, evtData);
  if (!_daughter2IsStable) _absDecDaughter2->fillWignerDs(fsMap, evtData);
  
  Vector4<double> all4Vec(0.,0.,0.,0.);
  Vector4<double> mother4Vec(0.,0.,0.,0.);
  Vector4<double> daughter2_4Vec(0.,0.,0.,0.);

  //fill all4Vec
  std::map<std::string, Vector4<double> >::iterator itMap;
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

  Vector4<double> daughter2HelMother(0.,0.,0.,0.);
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
  if (!_hasMotherPart && spinMother>1) lamMotherMax=1; //attention: this is only valid for pbar p or e+ e- reactions; must be moved to individual specific classes

 
  for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
    for (Spin lam12=-lam12Max; lam12<=lam12Max; ++lam12){
      double thePhi=0.;
      if(_hasMotherPart) thePhi=daughter2HelMother.Phi();
      evtData->WignerDsString[_wignerDKey][spinMother][lamMother][lam12]=Wigner_D(thePhi,daughter2HelMother.Theta(),0,spinMother,lamMother,lam12);
    }
  }
}

void AbsDecay::print(std::ostream& os) const{

  if(!_daughter1IsStable){
    os << "with further decay:";
    _absDecDaughter1->print(os);
  }
  if(!_daughter2IsStable){
    os << "with further decay:";
    _absDecDaughter2->print(os);
  }  

  os << "\n";
}
