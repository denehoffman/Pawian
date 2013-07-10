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
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/AbsDynamics.hh"


AbsDecay::AbsDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv) :
  _mother(mother)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(true)
  ,_motherIGJPCPtr(getIGJPCPtr(mother))
  ,_daughter1IGJPCPtr(getIGJPCPtr(daughter1))
  ,_daughter2IGJPCPtr(getIGJPCPtr(daughter2))
  ,_isospinClebschG(1.)
  ,_name(mother->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
  ,_massParamKey(_mother->name())
  ,_dynType("WoDynamics")
  ,_preFactor(1.)
  //  ,_dynKey(mother->name())
  ,_decPair1stChannel(make_pair(daughter1, daughter2))
  ,_env(theEnv)
  ,_gParity(mother->theGParity())
  ,_useIsospin(true)
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
  //  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+_motherJPCPtr->name()+FunctionUtils::particleListName(_finalStateParticles);
 _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+FunctionUtils::particleListName(_finalStateParticles);

 _idaughter1=Spin(_daughter1->twoIso(), 2);
 _i3daughter1=Spin(_daughter1->twoIso3(), 2);
 _idaughter2=Spin(_daughter2->twoIso(), 2);
 _i3daughter2=Spin(_daughter2->twoIso3(), 2);
 Spin Imother(_mother->twoIso(), 2);
 Spin I3mother(_mother->twoIso3(), 2);

 _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, _i3daughter2, Imother, I3mother);

 if(fabs(_isospinClebschG)<1.e-8){
   Warning << "no isospin coupling for decay " << _mother->name() << " to " << _daughter1->name() << " " << _daughter2->name() << endmsg;
   Warning << "Imother: " << Imother << "\tI3mother: " << I3mother << endmsg;
   Warning << "idaughter1: " << _idaughter1 << "\ti3daughter1: " << _i3daughter1 << endmsg;
   Warning << "idaughter2: " << _idaughter2 << "\ti3daughter2: " << _i3daughter2 << endmsg;
 }
}

AbsDecay::AbsDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName) :
  _mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(false)
  ,_motherIGJPCPtr(motherIGJPCPtr)
  ,_daughter1IGJPCPtr(getIGJPCPtr(daughter1))
  ,_daughter2IGJPCPtr(getIGJPCPtr(daughter2))
  ,_isospinClebschG(1.)
  ,_name(_motherIGJPCPtr->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_motherIGJPCPtr->jpcname()+"To"+daughter1->name()+"_"+daughter2->name())
  // ,_massParamKey(motherIGJPCPtr->name())
  ,_massParamKey(motherIGJPCPtr->jpcname())
  ,_dynType("WoDynamics")
  ,_preFactor(1.)
  //  ,_dynKey(motherJPCPtr->name())
  ,_decPair1stChannel(make_pair(daughter1, daughter2))
  ,_env(theEnv)
  ,_gParity(motherIGJPCPtr->G)
  ,_useIsospin(true)
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

  //  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+motherName;
  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+FunctionUtils::particleListName(_finalStateParticles);

  _daughter1->print(std::cout);
  _daughter2->print(std::cout);

  _idaughter1=Spin(_daughter1->twoIso(), 2);

  _i3daughter1=Spin(_daughter1->twoIso3(), 2);
  _idaughter2=Spin(_daughter2->twoIso(), 2);
  _i3daughter2=Spin(_daughter2->twoIso3(), 2);

  _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, _i3daughter2, motherIGJPCPtr->I, 0); //attention
}

AbsDecay::~AbsDecay(){
}


void AbsDecay::enableDynamics(std::string& dynString, std::vector<std::string>& additionalStringVec) {
  _dynType=dynString;

  if(additionalStringVec.size()==2){ //fill second decay channel (Flatte)
    Particle* firstParticle=_env->particleTable()->particle(additionalStringVec[0]);
    if(0==firstParticle){
      Alert << "particle with name\t" << additionalStringVec[0] << "\tnot available in the particle table" << endmsg;
      exit(1);
    }
    Particle* secondParticle=_env->particleTable()->particle(additionalStringVec[1]);
    if(0==secondParticle){
      Alert << "particle with name\t" << additionalStringVec[1] << "\tnot available in the particle table" << endmsg;
      exit(1);
    }
    _decPair2ndChannel=make_pair(firstParticle,secondParticle);
  }

  _absDynPtr=DynRegistry::instance()->getDynamics(shared_from_this());
}

void AbsDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, EvtData* evtData){
  int evtNo=evtData->evtNo;
  std::map<int, bool>::const_iterator it = _alreadyFilledMap.find(evtNo);
  if(it!=_alreadyFilledMap.end() &&  it->second) return; //already filled 

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

  Spin spinMother=_motherIGJPCPtr->J;
  Spin spinDaughter1=_daughter1IGJPCPtr->J;
  Spin spinDaughter2=_daughter2IGJPCPtr->J;
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

   _alreadyFilledMap[evtNo]=true;
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
