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
#include <sstream>
#include <algorithm>

#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDecayList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "ErrLogger/ErrLogger.hh"

AbsDecay::AbsDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelId, std::string typeName) :
  _typeName(typeName)
  ,_channelId(channelId)
  ,_mother(mother)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(true)
  ,_motherJPCPtr(getJPCPtr(mother))
  ,_motherIGJPCPtr(getIGJPCPtr(mother))
  ,_daughter1IGJPCPtr(getIGJPCPtr(daughter1))
  ,_daughter2IGJPCPtr(getIGJPCPtr(daughter2))
  ,_isospinClebschG(1.)
   ,_qR(BarrierFactor::qRDefault)
  ,_name(mother->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_name)
  ,_massParamKey(_mother->name())
   ,_prodParamKey("")
  ,_wignerDRefKey("default")
  ,_dynType("WoDynamics")
   ,_dynEnabled(false)
  ,_preFactor(1.)
   ,_pathParserFile("")
  //  ,_dynKey(mother->name())
  ,_decPair1stChannel(make_pair(daughter1, daughter2))
   ,_isDaughter1Photon(false)
   ,_isDaughter2Photon(false)
  ,_gParity(mother->theGParity())
  ,_useIsospin(true)
   ,_isProdAmp(false)
   ,_useProdBarrier(false)
   ,_massSumFsParticles(0.)
   ,_Lmin(0)
  ,_decLevel(decayLevel::noLevel)
  ,_prodChannelInfo(std::shared_ptr<ProdChannelInfo>(new ProdChannelInfo()))
{
  _absDecDaughter1=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter1);
  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    _finalStateParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter1.begin(), _finalStateParticlesDaughter1.end());
  }
  else{
    _finalStateParticles.push_back(daughter1);
    _finalStateParticlesDaughter1.push_back(daughter1);
    addToFsParticleNameList(daughter1->name());
  }

  _absDecDaughter2=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter2);

  if(0 != _absDecDaughter2){
    _daughter2IsStable=false;
    _finalStateParticlesDaughter2=_absDecDaughter2->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter2.begin(), _finalStateParticlesDaughter2.end());
  }
  else{
    _finalStateParticles.push_back(daughter2);
    _finalStateParticlesDaughter2.push_back(daughter2);
    addToFsParticleNameList(daughter2->name());
  }

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);
  //  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+_motherJPCPtr->name()+FunctionUtils::particleListName(_finalStateParticles);

 _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+FunctionUtils::particleListName(_finalStateParticles);
 _refKey=FunctionUtils::particleListName(_finalStateParticles);

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
 // if( (*daughter1) == *(GlobalEnv::instance()->particleTable()->particle("photon")) || (*daughter2) == *(GlobalEnv::instance()->particleTable()->particle("photon"))) disableIsospin();

 if(daughter1->twoJ() == 2 && daughter1->theParity() == -1 &&  daughter1->theCParity()==-1 && daughter1->mass() < 1.e-6){
   _isDaughter1Photon=true;   
 }
 if(daughter2->twoJ() == 2 && daughter2->theParity() == -1 &&  daughter2->theCParity()==-1 && daughter1->mass() < 1.e-6){
   _isDaughter2Photon=true;   
 }

 if( _isDaughter1Photon || _isDaughter2Photon) disableIsospin();
  
 //fill mass sum of final state particles
 std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    _massSumFsParticles+=(*itP)->mass();
  }
}

AbsDecay::AbsDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, std::string motherName, ChannelID channelId, std::string typeName) :
  _typeName(typeName)
  ,_channelId(channelId)
  ,_mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(false)
  ,_motherJPCPtr(std::shared_ptr<const jpcRes>(new jpcRes(motherIGJPCPtr->J, motherIGJPCPtr->P, motherIGJPCPtr->C)) )
  ,_motherIGJPCPtr(motherIGJPCPtr)
  ,_daughter1IGJPCPtr(getIGJPCPtr(daughter1))
  ,_daughter2IGJPCPtr(getIGJPCPtr(daughter2))
  ,_isospinClebschG(1.)
  ,_qR(BarrierFactor::qRDefault)
  ,_name(_motherIGJPCPtr->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_fitParamSuffix(_motherIGJPCPtr->jpcname()+"To"+daughter1->name()+"_"+daughter2->name())
  // ,_massParamKey(motherIGJPCPtr->name())
  ,_massParamKey(motherIGJPCPtr->jpcname())
  ,_prodParamKey(GlobalEnv::instance()->Channel(channelId)->channelTypeName()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_wignerDRefKey("default")
  ,_dynType("WoDynamics")
  ,_dynEnabled(false)
  ,_preFactor(1.)
  ,_pathParserFile("")
  //  ,_dynKey(motherJPCPtr->name())
  ,_decPair1stChannel(make_pair(daughter1, daughter2))
  ,_isDaughter1Photon(false)
  ,_isDaughter2Photon(false)
  ,_gParity(motherIGJPCPtr->G)
  ,_useIsospin(true)
  ,_isProdAmp(false)
  ,_useProdBarrier(false)
  ,_massSumFsParticles(0.)
  ,_Lmin(0)
  ,_decLevel(decayLevel::noLevel)
  ,_prodChannelInfo(std::shared_ptr<ProdChannelInfo>(new ProdChannelInfo()))
{
  _absDecDaughter1=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter1);

  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    _finalStateParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter1.begin(), _finalStateParticlesDaughter1.end());
  }
  else{
    _finalStateParticles.push_back(daughter1);
    _finalStateParticlesDaughter1.push_back(daughter1); 
    addToFsParticleNameList(daughter1->name());
  }

  _absDecDaughter2=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter2);

  if(0 != _absDecDaughter2){
    _daughter2IsStable=false;
    _finalStateParticlesDaughter2=_absDecDaughter2->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter2.begin(), _finalStateParticlesDaughter2.end());
  }
  else{
    _finalStateParticles.push_back(daughter2);
    _finalStateParticlesDaughter2.push_back(daughter2); 
    addToFsParticleNameList(daughter2->name());
  }

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  //  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+motherName;
  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+FunctionUtils::particleListName(_finalStateParticles);
  _refKey=FunctionUtils::particleListName(_finalStateParticles);
 
  _daughter1->print(std::cout);
  _daughter2->print(std::cout);

  _idaughter1=Spin(_daughter1->twoIso(), 2);

  _i3daughter1=Spin(_daughter1->twoIso3(), 2);
  _idaughter2=Spin(_daughter2->twoIso(), 2);
  _i3daughter2=Spin(_daughter2->twoIso3(), 2);

  if ( GlobalEnv::instance()->Channel(channelId)->channelType()==AbsChannelEnv::CHANNEL_PBARP 
       ||  GlobalEnv::instance()->Channel(channelId)->channelType()==AbsChannelEnv::CHANNEL_EPEM ){
    _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, _i3daughter2, motherIGJPCPtr->I, 0); //attention
  }
  else if (GlobalEnv::instance()->Channel(channelId)->channelType()==AbsChannelEnv::CHANNEL_GAMMAP){
    _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, _i3daughter2, motherIGJPCPtr->I, 1./2); //p gamma
  }
  else{
    Alert << "AbsDecay cannot be set up for channel type " <<  GlobalEnv::instance()->Channel(channelId)->channelType() << endmsg;
    exit(0);
  }

  //  if( (*daughter1) == *(GlobalEnv::instance()->particleTable()->particle("photon")) || (*daughter2) == *(GlobalEnv::instance()->particleTable()->particle("photon"))) disableIsospin();

 if(daughter1->twoJ() == 2 && daughter1->theParity() == -1 &&  daughter1->theCParity()==-1 && daughter1->mass() < 1.e-6){
   _isDaughter1Photon=true;   
 }
 if(daughter2->twoJ() == 2 && daughter2->theParity() == -1 &&  daughter2->theCParity()==-1 && daughter2->mass() < 1.e-6){
   _isDaughter2Photon=true;   
 }

 if( _isDaughter1Photon || _isDaughter2Photon) disableIsospin();

 //fill mass sum of final state particles
 std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    _massSumFsParticles+=(*itP)->mass();
  }
}

AbsDecay::~AbsDecay(){
}


void AbsDecay::enableDynamics(std::string& dynString, std::vector<std::string>& additionalStringVec) {
  if(_dynEnabled){
    Alert << "dynamics already enabled for " << name() << endmsg;
    exit(0);
  }
  _dynType=dynString;

  if(_dynType=="KMatrix") _pathParserFile=additionalStringVec[0];
  else if(_dynType=="Flatte"){ //fill second decay channel (Flatte)
    Particle* firstParticle=GlobalEnv::instance()->particleTable()->particle(additionalStringVec[0]);
    if(0==firstParticle){
      Alert << "particle with name\t" << additionalStringVec[0] << "\tnot available in the particle table" << endmsg;
      exit(1);
    }
    Particle* secondParticle=GlobalEnv::instance()->particleTable()->particle(additionalStringVec[1]);
    if(0==secondParticle){
      Alert << "particle with name\t" << additionalStringVec[1] << "\tnot available in the particle table" << endmsg;
      exit(1);
    }
    _decPair2ndChannel=make_pair(firstParticle,secondParticle);
  }
  else if(_dynType=="BlattWBarrier" || _dynType=="BlattWBarrierTensor" || _dynType=="BreitWignerBlattWRel" || _dynType=="BreitWignerBlattWTensorRel"){
    if(additionalStringVec.size()>0){

      _qR=stof(additionalStringVec[0]);
      if ( _qR <1.e-5 || _qR > 20.){
	Alert << "radius for barrier factor too high or too low qr=" << _qR << endmsg;
	exit(0);
      }
    }
    Info << "AmpName: " << name() << "  radius for barrier factor qr= " << _qR << endmsg;
  }
  
  _absDynPtr=DynRegistry::instance()->getDynamics(shared_from_this());
  _dynEnabled=true;
}

void AbsDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData){
  //  int evtNo=evtData->evtNo;
  //  std::map<int, bool>::const_iterator it = _alreadyFilledMap.find(evtNo);
  //  if(it!=_alreadyFilledMap.end() &&  it->second) return; //already filled

  std::vector<Particle*>::iterator itP;
  std::map<std::string, Vector4<double> >::iterator itMap;

  Vector4<double> all4Vec(0.,0.,0.,0.);
  Vector4<double> mother4Vec(0.,0.,0.,0.);
  Vector4<double> daughter1_4Vec(0.,0.,0.,0.);
  Vector4<double> daughter2_4Vec(0.,0.,0.,0.);

  //fill all4Vec
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill mother4Vec
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    mother4Vec+=itMap->second;
  }

  //fill daughter1
  for(itP=_finalStateParticlesDaughter1.begin(); itP!=_finalStateParticlesDaughter1.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    daughter1_4Vec+=itMap->second;
  }

  //fill daughter2
  for(itP=_finalStateParticlesDaughter2.begin(); itP!=_finalStateParticlesDaughter2.end(); ++itP){
    itMap=fsMap.find((*itP)->name());
    daughter2_4Vec+=itMap->second;
  }

  Vector4<double> daughter2HelMother(0.,0.,0.,0.);
  Vector4<double> daughter1HelMother(0.,0.,0.,0.);
  if(_hasMotherPart){
    if(fabs(mother4Vec==all4Vec)){
      daughter2HelMother=daughter2_4Vec;
      daughter2HelMother.Boost(daughter2HelMother); //is this correct????
      daughter1HelMother=daughter1_4Vec;
      daughter1HelMother.Boost(daughter1HelMother); //is this correct????
    }
    //    else daughter2HelMother=helicityVec(all4Vec, mother4Vec, daughter2_4Vec);
    else{
      daughter2HelMother=helicityVec(prodParticle4Vec, mother4Vec, daughter2_4Vec);
      daughter1HelMother=helicityVec(prodParticle4Vec, mother4Vec, daughter1_4Vec);
    }
  }
  else{
    daughter2HelMother=daughter2_4Vec;
    daughter2HelMother.Boost(mother4Vec);
    daughter1HelMother=daughter1_4Vec;
    daughter1HelMother.Boost(mother4Vec);
  }

  Spin spinMother=_motherIGJPCPtr->J;
  Spin spinDaughter1=_daughter1IGJPCPtr->J;
  Spin spinDaughter2=_daughter2IGJPCPtr->J;
  Spin lam12Max=spinDaughter1+spinDaughter2;
  if(lam12Max>spinMother) lam12Max=spinMother;

  Spin lamMotherMax=spinMother;
  if (!_hasMotherPart && spinMother>1){
    if ( GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_PBARP 
	 ||  GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_EPEM ){
      lamMotherMax=1; //attention
    }
    else if (GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_GAMMAP){
      lamMotherMax=3./2;
    }
  }

  for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
    for (Spin lam12=-lam12Max; lam12<=lam12Max; ++lam12){
      double thePhi=0.;
      if(_hasMotherPart) thePhi=daughter2HelMother.Phi();
      Id3StringType IdSpinMotherLamMotherLam12=FunctionUtils::spin3Index(spinMother, lamMother, lam12);
      
      std::map<Id3StringType, complex<double> >::iterator found = evtData->WignerDStringStringId[_wignerDKey][_wignerDRefKey].find(IdSpinMotherLamMotherLam12);
      if(found != evtData->WignerDStringStringId[_wignerDKey][_wignerDRefKey].end()){
	continue;
      }
      
      if (GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_EPEM && whichDecayLevel()==decayLevel::isProdAmp){
	evtData->WignerDStringStringId[_wignerDKey][_wignerDRefKey][IdSpinMotherLamMotherLam12]=Wigner_D(thePhi,daughter1HelMother.Theta(),0,spinMother,lamMother,lam12);
      } 
      else evtData->WignerDStringStringId[_wignerDKey][_wignerDRefKey][IdSpinMotherLamMotherLam12]=Wigner_D(thePhi,daughter2HelMother.Theta(),0,spinMother,lamMother,lam12);
      if(evtData->WignerDStringStringId[_wignerDKey][_wignerDRefKey][IdSpinMotherLamMotherLam12] != evtData->WignerDStringStringId[_wignerDKey][_wignerDRefKey][IdSpinMotherLamMotherLam12]){
	DebugMsg << "WignerDsString nan:\t" << evtData->WignerDsString[_wignerDKey][spinMother][lamMother][lam12] << endmsg;
	DebugMsg << "daughter2HelMother.Theta():\t" << daughter2HelMother.Theta() << endmsg;
	DebugMsg << "thePhi:\t" << thePhi << endmsg;
	DebugMsg << "daughter2_4Vec:\t" << daughter2_4Vec << endmsg;
	DebugMsg << "daughter2HelMother:\t" << daughter2HelMother << endmsg;
	DebugMsg << "mother4Vec:\t" << mother4Vec << "\tP:" << mother4Vec.P() << endmsg;
      }
    }
  }

  if (!_daughter1IsStable){
    _absDecDaughter1->fillWignerDs(fsMap, mother4Vec, evtData); 
  }
  if (!_daughter2IsStable){
    _absDecDaughter2->fillWignerDs(fsMap, mother4Vec, evtData);
  }

  bool fillqVals=false;
  if(_isProdAmp && _useProdBarrier) fillqVals=true;
  else if(0!=_absDynPtr){
    if(_absDynPtr->type()=="BlattWBarrierDynamics" || _absDynPtr->type()=="BlattWBarrierTensorDynamics") fillqVals=true; 
  }

  if(fillqVals){
       double qVal=daughter2HelMother.P();
    double qValNorm=breakupMomQ(mother4Vec.M(), massSumFsParticlesDec1(), massSumFsParticlesDec2()).real();
    evtData->DoubleString[_wignerDKey]=qVal;
    evtData->DoubleString[_wignerDKey+"qNorm"] = qValNorm;
  } 
  //   _alreadyFilledMap[evtNo]=true;
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


 void AbsDecay::resetFilledMap() {
   _alreadyFilledMap.clear();

   if (!_daughter1IsStable){
     _absDecDaughter1->resetFilledMap();
   }
   if (!_daughter2IsStable){
     _absDecDaughter2->resetFilledMap();
   }
   
}

void AbsDecay::enableProdBarrier(){
  if(_dynEnabled){
    Alert << "dynamics already enabled for " << name() << endmsg;
    exit(1);
  }
  if(!_prodChannelInfo->isProductionChannel()){
    Warning << name() << " is not a production amplitide! Barrier factors for the production can not be enabled!" << endmsg;
    return;
  }
  if(!_prodChannelInfo->withProdBarrier()){
    Warning << name() << "production barrier disabled" << endmsg;
    return;
  }
  _useProdBarrier=true;
  if(_prodChannelInfo->prodBarrierType()=="BlattWBarrier"){
    _dynType="BlattWBarrier";
    _qR=_prodChannelInfo->qRPod();
  }
  else if(_prodChannelInfo->prodBarrierType()=="RadM1"){
    _dynType="RadM1";
  }

  Info << "Barrier factors for production amplitude " << name() << " enabled!" << endmsg;
  _absDynPtr=DynRegistry::instance()->getDynamics(shared_from_this()); 
  _dynEnabled=true;
}

double AbsDecay::massSumFsParticlesDec1(){
  if (!_daughter1IsStable) return _absDecDaughter1->massSumFsParticles();
  else return _daughter1->mass();
}

double AbsDecay::massSumFsParticlesDec2(){
  if (!_daughter2IsStable) return _absDecDaughter2->massSumFsParticles();
  return _daughter2->mass();
}

void AbsDecay::extractLmin(){
  std::vector< std::shared_ptr<const LScomb> > LSDecAmps;
  if (_useIsospin){
    Spin currentGParity=_motherIGJPCPtr->G;
    int daughter1GParity=_daughter1->theGParity();
    int daughter2GParity=_daughter2->theGParity();

    if( fabs(currentGParity)==1 && fabs(daughter1GParity)==1 && fabs(daughter2GParity)==1){
      validLS( _motherIGJPCPtr, _daughter1, _daughter2, LSDecAmps, true, _gParity, true);
    }
    else{
      validLS( _motherIGJPCPtr, _daughter1, _daughter2, LSDecAmps);
    }
  }
  else{
    validLS( _motherIGJPCPtr, _daughter1, _daughter2, LSDecAmps);
  }
 
  _Lmin=100;
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=LSDecAmps.begin(); it!=LSDecAmps.end(); ++it){
    int currentL=(*it)->L;
    if(_Lmin>currentL) _Lmin=currentL;
  }
   
}

void AbsDecay::setDecayLevel(decLevel theLevel){
  if (_decLevel==decayLevel::noLevel){
    _decLevel=theLevel;
  }
  else if(_decLevel != theLevel) _decLevel=decayLevel::severalLevels;
  else return;
  Info << name() << " set decay level to " << _decLevel << endmsg; 
}

void AbsDecay::setDecayLevelTree(decLevel theLevel, std::shared_ptr<AbsDecay> motherDecPtr, std::shared_ptr<AbsDecay> prodDecPtr){
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

  if (!_daughter1IsStable){
    _absDecDaughter1->setDecayLevelTree(decayLevel(theLevel+1), shared_from_this(), prodDecPtr);
  }
  else{
    prodDecPtr->addToFsParticleNameList(_daughter1->name());
    motherDecPtr->addToFsParticleNameList(_daughter1->name());
  }

  if (!_daughter2IsStable){
    _absDecDaughter2->setDecayLevelTree(decayLevel(theLevel+1), shared_from_this(), prodDecPtr);
  }
  else{
    prodDecPtr->addToFsParticleNameList(_daughter2->name());
    motherDecPtr->addToFsParticleNameList(_daughter2->name());
  }
}

void AbsDecay::addToFsParticleNameList(const std::string& name){
  std::vector<std::string>::iterator it;
  bool alreadyThere=false;
  for(it = _fsParticleNameList.begin(); it!=_fsParticleNameList.end(); ++it){
    if (name==(*it)) alreadyThere=true;
  }
  if( !alreadyThere ) _fsParticleNameList.push_back(name);
  
}

void  AbsDecay::setWigDRefKey(std::string& ref){
  if(ref == _wignerDRefKey) return;
  if (_wignerDRefKey != "default"){
    Alert << "two different reference keys for the WignerD functions not allowed!!!"
	  << "\n which are: " <<  _wignerDRefKey << " and: " << ref 
	  << "\n clone particle(s) for at least one decay chain" << endmsg;
    exit(1);
  }
  _wignerDRefKey=ref;

  Info << _name << ":\twignerDKey= " << _wignerDKey << "\twignerDRefKey= " << _wignerDRefKey << endmsg;

  if (!_daughter1IsStable){
    _absDecDaughter1->setWigDRefKey(_refKey);
  }
  if (!_daughter2IsStable){
    _absDecDaughter2->setWigDRefKey(_refKey);
  }
}
