//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universit??t Bochum 				  //
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
#include "qft++Extension/PawianUtils.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/PawianConstants.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/ProdChannelInfo.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "ErrLogger/ErrLogger.hh"
#include "ConfigParser/ParserBase.hh"

AbsDecay::AbsDecay(Particle* mother, Particle* daughter1, Particle* daughter2, 
		   ChannelID channelId, std::string typeName) :
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
  ,_nameId(0)
  ,_fitParamSuffix(_name)
  ,_massParamKey(_mother->name())
  ,_massParamId(IdStringMapRegistry::instance()->keyStringId("grandMaAndMassParKey", _massParamKey))
   ,_prodParamKey("")
  ,_wignerDId(0)
  ,_wignerDqNormId(0)
  ,_wignerDRefKey("default")
  ,_wigDWigDRefId(0)
  ,_dynType("WoDynamics")
   ,_dynEnabled(false)
  ,_preFactor(1.)
   ,_pathParserFile("")
  ,_projectionParticleNames("")
  ,_decPair1stChannel(make_pair(daughter1, daughter2))
   ,_isDaughter1Photon(false)
   ,_isDaughter2Photon(false)
  ,_gParity(mother->theGParity())
  ,_useIsospin(true)
  ,_isWeakDecay(false)
   ,_isProdAmp(false)
   ,_useProdBarrier(false)
   ,_massSumFsParticles(0.)
   ,_Lmin(0)
  ,_decLevel(decayLevel::noLevel)
  ,_motherName("")
  ,_prodChannelInfo(std::shared_ptr<ProdChannelInfo>(new ProdChannelInfo()))
  ,_useMultipoleGeneral(GlobalEnv::instance()->Channel(_channelId)->parser()->useMultipoleGeneral())
{
  _absDecDaughter1=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter1);
  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    _finalStateParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter1.begin(), 
				_finalStateParticlesDaughter1.end());
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
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter2.begin(), 
				_finalStateParticlesDaughter2.end());
  }
  else{
    _finalStateParticles.push_back(daughter2);
    _finalStateParticlesDaughter2.push_back(daughter2);
    addToFsParticleNameList(daughter2->name());
  }
  
  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2) +
    "_"+FunctionUtils::particleListName(_finalStateParticles);
  _refKey=FunctionUtils::particleListName(_finalStateParticles);
  
  _idaughter1=Spin(_daughter1->twoIso(), 2);
  _i3daughter1=Spin(_daughter1->twoIso3(), 2);
  _idaughter2=Spin(_daughter2->twoIso(), 2);
  _i3daughter2=Spin(_daughter2->twoIso3(), 2);
  
  //check whether it's a weak decay 
  if( _mother->strange() != (_daughter1->strange()+_daughter2->strange())) _isWeakDecay=true;
  else if ( _mother->charm() != (_daughter1->charm()+_daughter2->charm())) _isWeakDecay=true;
  if(_isWeakDecay) disableIsospin(); 
  else {  
    Spin Imother(_mother->twoIso(), 2);
    Spin I3mother(_mother->twoIso3(), 2);
    
    _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, _i3daughter2, Imother, I3mother);

    //check z-component of the isospin
    if( I3mother != (_i3daughter1+_i3daughter2)){
      Alert << "electric charge is not conserved for the decay " << _mother->name() << " to " 
	    << _daughter1->name() << " " << _daughter2->name() << endmsg;
      Alert << "I3(mother): " << I3mother << " != " << "I3(daughter1): " <<  _i3daughter1 
	    << " + " << "I3(daughter2): " <<  _i3daughter2 << endmsg;
      exit(1); 
    } 
  
    if(fabs(_isospinClebschG)<1.e-8){
      WarningMsg << "no isospin coupling for decay " << _mother->name() 
		 << " to " << _daughter1->name() << " " << _daughter2->name() << endmsg;
      WarningMsg << "Imother: " << Imother << "\tI3mother: " << I3mother << endmsg;
      WarningMsg << "idaughter1: " << _idaughter1 << "\ti3daughter1: " << _i3daughter1 << endmsg;
      WarningMsg << "idaughter2: " << _idaughter2 << "\ti3daughter2: " << _i3daughter2 << endmsg;
    }
    
    if(daughter1->twoJ() == 2 && daughter1->theParity() == -1 &&  daughter1->theCParity()==-1 
       && daughter1->mass() < 1.e-6){
      _isDaughter1Photon=true;   
    }
    if(daughter2->twoJ() == 2 && daughter2->theParity() == -1 &&  daughter2->theCParity()==-1 
       && daughter1->mass() < 1.e-6){
      _isDaughter2Photon=true;   
    }
    
    if( _isDaughter1Photon || _isDaughter2Photon) disableIsospin();
  }

  //fill mass sum of final state particles
  std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    _massSumFsParticles+=(*itP)->mass();
  }

}

AbsDecay::AbsDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, 
		   std::string motherName, ChannelID channelId, std::string typeName) :
  _typeName(typeName)
  ,_channelId(channelId)
  ,_mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(0)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(false)
  ,_motherJPCPtr(std::shared_ptr<const jpcRes>(new jpcRes(motherIGJPCPtr->J, 
							  motherIGJPCPtr->P, 
							  motherIGJPCPtr->C)) )
  ,_motherIGJPCPtr(motherIGJPCPtr)
  ,_daughter1IGJPCPtr(getIGJPCPtr(daughter1))
  ,_daughter2IGJPCPtr(0)
  ,_isospinClebschG(1.)
  ,_qR(BarrierFactor::qRDefault)
  ,_name(_motherIGJPCPtr->name()+"To"+daughter1->name())
  ,_nameId(0)
  ,_fitParamSuffix(_motherIGJPCPtr->jpcname()+"To"+daughter1->name())
  ,_massParamKey(motherIGJPCPtr->jpcname())
  ,_massParamId(IdStringMapRegistry::instance()->keyStringId("grandMaAndMassParKey", _massParamKey))
  ,_prodParamKey(GlobalEnv::instance()->Channel(channelId)->channelTypeName()+"To"+daughter1->name())
  ,_wignerDId(0)
  ,_wignerDqNormId(0)
  ,_wignerDRefKey("default")
  ,_wigDWigDRefId(0)
  ,_dynType("WoDynamics")
  ,_dynEnabled(false)
  ,_preFactor(1.)
  ,_pathParserFile("")
  ,_projectionParticleNames("")
  ,_isDaughter1Photon(false)
  ,_isDaughter2Photon(false)
  ,_gParity(motherIGJPCPtr->G)
  ,_useIsospin(true)
  ,_isWeakDecay(false)
  ,_isProdAmp(true)
  ,_useProdBarrier(false)
  ,_massSumFsParticles(0.)
  ,_Lmin(0)
  ,_decLevel(decayLevel::isProdAmp)
  ,_motherName(motherName)
  ,_prodChannelInfo(std::shared_ptr<ProdChannelInfo>(new ProdChannelInfo()))
  ,_useMultipoleGeneral(GlobalEnv::instance()->Channel(_channelId)->parser()->useMultipoleGeneral())
{
  _absDecDaughter1=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter1);

  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    _finalStateParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter1.begin(), 
				_finalStateParticlesDaughter1.end());
  }
  else if(_daughter1->name() == "GamGam"){
   _daughter1IsStable=true;
  }
  else{
    Alert << "formed particle " << _daughter1->name()  << " can not be a stable particle!!!" << endmsg;
    exit(1);
  }

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_" +
    FunctionUtils::particleListName(_finalStateParticles);
  _refKey=FunctionUtils::particleListName(_finalStateParticles);
 
  if (ErrLogger::instance().level()>logging::log_level::NOTICE) _daughter1->print(std::cout);

  _idaughter1=Spin(_daughter1->twoIso(), 2);

  _i3daughter1=Spin(_daughter1->twoIso3(), 2);

  //fill mass sum of final state particles
  std::vector<Particle*>::iterator itP;
  for(itP = _finalStateParticles.begin(); itP != _finalStateParticles.end(); ++itP){
    _massSumFsParticles+=(*itP)->mass();
  }
}

AbsDecay::AbsDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, 
		   Particle* daughter2, std::string motherName, ChannelID channelId, std::string typeName) :
  _typeName(typeName)
  ,_channelId(channelId)
  ,_mother(0)
  ,_daughter1(daughter1)
  ,_daughter2(daughter2)
  ,_daughter1IsStable(true)
  ,_daughter2IsStable(true)
  ,_hasMotherPart(false)
  ,_motherJPCPtr(std::shared_ptr<const jpcRes>(new jpcRes(motherIGJPCPtr->J, 
							  motherIGJPCPtr->P, 
							  motherIGJPCPtr->C)) )
  ,_motherIGJPCPtr(motherIGJPCPtr)
  ,_daughter1IGJPCPtr(getIGJPCPtr(daughter1))
  ,_daughter2IGJPCPtr(getIGJPCPtr(daughter2))
  ,_isospinClebschG(1.)
  ,_qR(BarrierFactor::qRDefault)
  ,_name(_motherIGJPCPtr->name()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_nameId(0)
  ,_fitParamSuffix(_motherIGJPCPtr->jpcname()+"To"+daughter1->name()+"_"+daughter2->name())
  ,_massParamKey(motherIGJPCPtr->jpcname())
  ,_massParamId(IdStringMapRegistry::instance()->keyStringId("grandMaAndMassParKey", _massParamKey))
  ,_prodParamKey(GlobalEnv::instance()->Channel(channelId)->channelTypeName()+"To" +
		 daughter1->name()+"_"+daughter2->name())
  ,_wignerDId(0)
  ,_wignerDqNormId(0)
  ,_wignerDRefKey("default")
  ,_wigDWigDRefId(0)
  ,_dynType("WoDynamics")
  ,_dynEnabled(false)
  ,_preFactor(1.)
  ,_pathParserFile("")
  ,_projectionParticleNames("")
  //  ,_dynKey(motherJPCPtr->name())
  ,_decPair1stChannel(make_pair(daughter1, daughter2))
  ,_isDaughter1Photon(false)
  ,_isDaughter2Photon(false)
  ,_gParity(motherIGJPCPtr->G)
  ,_useIsospin(true)
  ,_isWeakDecay(false)
  ,_isProdAmp(false)
  ,_useProdBarrier(false)
  ,_massSumFsParticles(0.)
  ,_Lmin(0)
  ,_decLevel(decayLevel::noLevel)
  ,_prodChannelInfo(std::shared_ptr<ProdChannelInfo>(new ProdChannelInfo()))
  ,_useMultipoleGeneral(GlobalEnv::instance()->Channel(_channelId)->parser()->useMultipoleGeneral())
{
  _absDecDaughter1=GlobalEnv::instance()->Channel(_channelId)->absDecayList()->decay(_daughter1);

  if(0 != _absDecDaughter1){
    _daughter1IsStable=false;
    _finalStateParticlesDaughter1=_absDecDaughter1->finalStateParticles();
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter1.begin(), 
				_finalStateParticlesDaughter1.end());
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
    _finalStateParticles.insert(_finalStateParticles.end(), _finalStateParticlesDaughter2.begin(), 
				_finalStateParticlesDaughter2.end());
  }
  else{
    _finalStateParticles.push_back(daughter2);
    _finalStateParticlesDaughter2.push_back(daughter2); 
    addToFsParticleNameList(daughter2->name());
  }

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  //  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_"+motherName;
  _wignerDKey=FunctionUtils::particleListName(_finalStateParticlesDaughter2)+"_" +
    FunctionUtils::particleListName(_finalStateParticles);
  _refKey=FunctionUtils::particleListName(_finalStateParticles);
 
  if (ErrLogger::instance().level()>logging::log_level::NOTICE) _daughter1->print(std::cout);
  if (ErrLogger::instance().level()>logging::log_level::NOTICE) _daughter2->print(std::cout);

  _idaughter1=Spin(_daughter1->twoIso(), 2);

  _i3daughter1=Spin(_daughter1->twoIso3(), 2);
  _idaughter2=Spin(_daughter2->twoIso(), 2);
  _i3daughter2=Spin(_daughter2->twoIso3(), 2);

  if ( GlobalEnv::instance()->Channel(channelId)->channelType()==AbsChannelEnv::CHANNEL_PBARP 
       ||  GlobalEnv::instance()->Channel(channelId)->channelType()==AbsChannelEnv::CHANNEL_EPEM ){
    _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, 
			     _i3daughter2, motherIGJPCPtr->I, 0); //attention
    if ( (_i3daughter1+_i3daughter2) != Spin(0) ){
      Alert << "electric charge is not conserved for the production of the two daughters " 
	    << _daughter1->name() << " and " << _daughter2->name() << endmsg;
      Alert << "I3(daughter1): " <<  _i3daughter1 << " + " << "I3(daughter2): " 
	    <<  _i3daughter2 << "  !=  0" << endmsg;
      exit(1);
    }
  }
  else if (GlobalEnv::instance()->Channel(channelId)->channelType()==AbsChannelEnv::CHANNEL_GAMMAP){
    _isospinClebschG=Clebsch(_idaughter1, _i3daughter1, _idaughter2, 
			     _i3daughter2, motherIGJPCPtr->I, 1./2); //p gamma
    if ( (_i3daughter1+_i3daughter2) != Spin(1./2) ){
      Alert << "electric charge is not conserved for the production of the two daughters " 
	    << _daughter1->name() << " and " << _daughter2->name() << endmsg;
      Alert << "I3(daughter1): " <<  _i3daughter1 << " + " << "I3(daughter2): " <<  
	_i3daughter2 << "  !=  1/2" << endmsg;
      exit(1);
    }
  }
  else{
    Alert << "AbsDecay cannot be set up for channel type " 
	  <<  GlobalEnv::instance()->Channel(channelId)->channelType() << endmsg;
    exit(0);
  }

 if(daughter1->twoJ() == 2 && daughter1->theParity() == -1 &&  daughter1->theCParity()==-1 
    && daughter1->mass() < 1.e-6){
   _isDaughter1Photon=true;   
 }
 if(daughter2->twoJ() == 2 && daughter2->theParity() == -1 &&  daughter2->theCParity()==-1 
    && daughter2->mass() < 1.e-6){
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

  if(_dynType=="KMatrix" || _dynType=="TMatrix" || _dynType=="TMatrixCompare" || _dynType=="FVectorCompare" 
     || _dynType=="FVectorIntensity" || _dynType=="Omnes"){
    _pathParserFile=additionalStringVec[0];
    if (additionalStringVec.size() == 3) 
      _projectionParticleNames=additionalStringVec.at(1)+"\t"+additionalStringVec.at(2);
  }
  else if(_dynType=="Flatte"){ //fill second decay channel (Flatte)
    Particle* firstParticle=GlobalEnv::instance()->particleTable()->particle(additionalStringVec[0]);
    if(0==firstParticle){
      Alert << "particle with name\t" << additionalStringVec[0] 
	    << "\tnot available in the particle table" << endmsg;
      exit(1);
    }
    Particle* secondParticle=GlobalEnv::instance()->particleTable()->particle(additionalStringVec[1]);
    if(0==secondParticle){
      Alert << "particle with name\t" << additionalStringVec[1] 
	    << "\tnot available in the particle table" << endmsg;
      exit(1);
    }
    _decPair2ndChannel=make_pair(firstParticle,secondParticle);
  }
  else if(_dynType=="BlattWBarrier" || _dynType=="BlattWBarrierTensor" 
	  || _dynType=="BreitWignerBlattWRel" || _dynType=="BreitWignerBlattWTensorRel"){
    if(additionalStringVec.size()>0){

      _qR=stof(additionalStringVec[0]);
      if ( _qR <1.e-5 || _qR > 20.){
	Alert << "radius for barrier factor too high or too low qr=" << _qR << endmsg;
	exit(0);
      }
    }
    InfoMsg << "AmpName: " << name() << "  radius for barrier factor qr= " << _qR << endmsg;
  }
  else if(_dynType=="LinearDynamics"){
    if (additionalStringVec.size() != 1){
      Alert << "for dynamics LinearDynamics value for reference mass must be set" << endmsg;
	exit(0);
    }
    _refMassLinearDyn=stof(additionalStringVec.at(0));
    InfoMsg << "AmpName: " << name() << "  reference mass for linear dynamics = " 
	    << _refMassLinearDyn << endmsg;
  }

  _absDynPtr=DynRegistry::instance()->getDynamics(shared_from_this());
  _dynEnabled=true;
}

void AbsDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, 
			    Vector4<double>& prodParticle4Vec, EvtData* evtData){
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
  Vector4<double> motherRefVec;


  if(_hasMotherPart){
    if(whichDecayLevel()==decayLevel::isProdAmp){
      Alert << "this amplitude " << name() <<" has got a mother and is a production amplitude!!!" << endmsg;
      exit(1); 
    }
    if(whichDecayLevel()==decayLevel::firstLevel){
      motherRefVec=Vector4<double>(GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().E(),
				   GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().Px(),
				   GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().Py(),
				   GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().Pz());
;
 if( GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_GG ){
   // real photons: quantization axis is assumed to be the z-axis 
	motherRefVec=Vector4<double>(10., 3., 0., 0.);
	prodParticle4Vec=Vector4<double>(GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().E(),
				   0.,
				   0.,
				   GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().P());	
      }
    }
    else if(whichDecayLevel()==decayLevel::secondLevel){
      motherRefVec=all4Vec;
      if( GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_GG ){
	// real photons: quantization axis is assumed to be the z-axis 
	motherRefVec=Vector4<double>(GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().E(),
				     0.,
				     0.,
				     GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().P());	

      }
    }
    else{
      std::vector< std::shared_ptr<AbsDecay> > motherAmpList=motherAmpRefList();
      if(motherAmpList.size()>1 || motherAmpList.size()==0){
        Alert << "decay level " << whichDecayLevel() << " for the decay " << name() << " contains 0 or more than 1 production amplitudes! Production plane cannot be set properly" << endmsg;
        Alert << "motherAmpList.size(): " << motherAmpList.size() << endmsg;
        exit(0);
      }
      
       std::vector< std::shared_ptr<AbsDecay> > motherMotherAmpList=motherAmpList.at(0)->motherAmpRefList();
      if(motherMotherAmpList.size()>1 || motherMotherAmpList.size()==0){
        Alert << "decay level " << whichDecayLevel() << " for the decay " << name() << " contains 0 or more than 1 production reference amplitudes! Production plane cannot be set properly" << endmsg;
        Alert << "motherMotherAmpList.size(): " << motherMotherAmpList.size() << endmsg;
        exit(0);
      }
 
      Vector4<double> motherRefDaughter1_4Vec;
      std::vector<Particle*> motherRefDaughter1List=motherMotherAmpList.at(0)->finalStateParticlesDaughter1();
      for(itP=motherRefDaughter1List.begin(); itP!=motherRefDaughter1List.end(); ++itP){
	itMap=fsMap.find((*itP)->name());
	motherRefDaughter1_4Vec+=itMap->second;
      }
      
      //fill daughter2
      Vector4<double> motherRefDaughter2_4Vec;
      std::vector<Particle*> motherRefDaughter2List=motherMotherAmpList.at(0)->finalStateParticlesDaughter2();
      for(itP=motherRefDaughter2List.begin(); itP!=motherRefDaughter2List.end(); ++itP){
	itMap=fsMap.find((*itP)->name());
	motherRefDaughter2_4Vec+=itMap->second;
      }
      
      motherRefVec=motherRefDaughter1_4Vec+motherRefDaughter2_4Vec;
      //InfoMsg << "motherRefVec: " << motherRefVec << endmsg;
      //InfoMsg << "prodParticle4Vec: " << prodParticle4Vec << endmsg;
    }  

    daughter2HelMother=KinUtils::heliVec(motherRefVec, prodParticle4Vec, mother4Vec, daughter2_4Vec);
    daughter1HelMother=KinUtils::heliVec(motherRefVec, prodParticle4Vec, mother4Vec, daughter1_4Vec);
  }
  else{
    if(whichDecayLevel()!=decayLevel::isProdAmp){
      Alert << "this amplitude " << name() 
	    << " hasn't got a mother and is not a production amplitude!!!" << endmsg;
      exit(1); 
    }
    daughter2HelMother=daughter2_4Vec;
    daughter1HelMother=daughter1_4Vec;
    if( fabs(mother4Vec.P()) > 1.e-6 ){
      Vector4<double> defaultMotherRefVec=Vector4<double>(10., 3., 0., 0.);
      Vector4<double> defaultRefVec =
	Vector4<double>(GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().E(),
			GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().Px(),
			GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().Py(),
			GlobalEnv::instance()->Channel(_channelId)->projectile4Vec().Pz());
      
      daughter2HelMother=KinUtils::heliVec(defaultMotherRefVec, defaultRefVec, mother4Vec, daughter2_4Vec);
      daughter1HelMother=KinUtils::heliVec(defaultMotherRefVec, defaultRefVec, mother4Vec, daughter1_4Vec);
    }
  }

  Spin spinMother=_motherIGJPCPtr->J;
  Spin spinDaughter1=_daughter1IGJPCPtr->J;
  Spin spinDaughter2=_daughter2IGJPCPtr->J;
  Spin lam12Max=spinDaughter1+spinDaughter2;
  if(lam12Max>spinMother) lam12Max=spinMother;

  Spin lamMotherMax=spinMother;
  if (whichDecayLevel()==decayLevel::isProdAmp && spinMother>1){
    if ( GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_PBARP 
	 ||  GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_EPEM ){
      lamMotherMax=1; //attention
    }
    else if (GlobalEnv::instance()->Channel(_channelId)->channelType()==AbsChannelEnv::CHANNEL_GAMMAP){
      lamMotherMax=3./2;
    }
  }

  double theTheta=daughter1HelMother.Theta();
  double thePhi=daughter1HelMother.Phi();
  
  for (Spin lamMother=-lamMotherMax; lamMother<=lamMotherMax; ++lamMother){
    for (Spin lam12=-lam12Max; lam12<=lam12Max; ++lam12){
      Id3StringType IdSpinMotherLamMotherLam12=FunctionUtils::spin3Index(spinMother, lamMother, lam12);
      std::map<Id3StringType, complex<double> >::iterator found = 
	evtData->WignerDIdId3[_wigDWigDRefId].find(IdSpinMotherLamMotherLam12);
      if(found != evtData->WignerDIdId3[_wigDWigDRefId].end()){
	continue;
      }

      if(type()=="IsobarHeliMultipoleDecay"){
        //order of the fs particles are fixed here. daughter1 is the photon
        //thePhi=PawianConstants::pi+thePhi;
        //theTheta=PawianConstants::pi-daughter1HelMother.Theta();
      	//evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12] =
        //  Wigner_D(thePhi,theTheta,0,spinMother,lamMother,lam12);

        evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12] =
          conj(Wigner_D(thePhi,theTheta,0,spinMother,lamMother,lam12));
      }
      // if(type()=="IsobarHeliMultipoleDecay"){
      // 	evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12] =
      // 	  conj(Wigner_D(thePhi,theTheta,0,spinMother,lamMother,lam12));
      // }
      
      else if (GlobalEnv::instance()->Channel(_channelId)->parser()->productionFormalism()=="HeliMultipole"){
	evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12] =
	  conj(Wigner_D(thePhi,theTheta,0,spinMother,lamMother,lam12));
      }
      else{
	evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12] =
	  Wigner_D(thePhi,theTheta,0,spinMother,lamMother,lam12);
      }

      
      if(evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12].real() != 
	 evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12].real()){
	Alert << "WignerD function of event No: " << evtData->evtNo << " is nan!!! " << endmsg;
	Alert << "Set it manually to complex<double>(1.e-10, 1.e-10)" << endmsg;
	evtData->WignerDIdId3[_wigDWigDRefId][IdSpinMotherLamMotherLam12]= complex<double>(1.e-10, 1.e-10);
	Alert << "daughter1HelMother phi: " << daughter1HelMother.Phi() << "\ttheta: " 
	      << daughter1HelMother.Theta() << endmsg;
	Alert << "spinMother: " << spinMother << "\tlam12: " << lam12 << endmsg;
	Alert << "daughter1HelMother: Mass: " << daughter1HelMother.Mass() 
	      << "\tE: " << daughter1HelMother.E()
	      << "\tPx: " << daughter1HelMother.Px()
	      << "\tPy: " << daughter1HelMother.Py()
	      << "\tPz: " << daughter1HelMother.Pz()  
	      << endmsg;
	
	Alert << "motherRefVec: " << motherRefVec << endmsg;
	Alert << "prodParticle4Vec: " << prodParticle4Vec << endmsg;
	Alert << "mother4Vec: " << mother4Vec << endmsg;
	Alert << "daughter1_4Vec: " << daughter1_4Vec << endmsg;
	Alert << "daughter2_4Vec: " << daughter2_4Vec << endmsg;
       	exit(1); 
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
    if(_absDynPtr->type()=="BlattWBarrierDynamics" || _absDynPtr->type()=="BlattWBarrierTensorDynamics") 
      fillqVals=true; 
  }

  if(fillqVals){
       double qVal=daughter2HelMother.P();
    double qValNorm=PawianQFT::breakupMomQDefault(mother4Vec.M(), massSumFsParticlesDec1(), 
						  massSumFsParticlesDec2()).real();
    evtData->DoubleMassId[_wignerDqId]=qVal;
    evtData->DoubleMassId[_wignerDqNormId] = qValNorm;
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
    WarningMsg << name() << " is not a production amplitide! Barrier factors for the production "
	       << "can not be enabled!" << endmsg;
    return;
  }
  if(!_prodChannelInfo->withProdBarrier()){
    WarningMsg << name() << "production barrier disabled" << endmsg;
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
  else if(_prodChannelInfo->prodBarrierType()=="RadM1KEDR"){
    _dynType="RadM1KEDR";
  }
  else if(_prodChannelInfo->prodBarrierType()=="SExpDynamics"){
    _dynType="SExpDynamics";
  }
  InfoMsg << "_prodChannelInfo->prodBarrierType(): " << _prodChannelInfo->prodBarrierType() << endmsg;
  InfoMsg << "Barrier factors for production amplitude " << name() << " enabled!" << endmsg;
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
    if(_decLevel==decayLevel::isProdAmp) _hasMotherPart=false;
  }
  else if(_decLevel != theLevel) _decLevel=decayLevel::severalLevels;
  else return;
  InfoMsg << name() << " set decay level to " << _decLevel << endmsg; 
}

void AbsDecay::setDecayLevelTree(decLevel theLevel, std::shared_ptr<AbsDecay> motherDecPtr, 
				 std::shared_ptr<AbsDecay> prodDecPtr){
  setDecayLevel(theLevel);
  if(prodDecPtr->whichDecayLevel() != AbsDecay::decayLevel::isProdAmp){
    Alert << "prodDecPtr with the name " << prodDecPtr->name() 
	  << " is not a production amplitude!!!" << endmsg;
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
    Alert << "for decay: " << name() <<" two different reference keys for the WignerD functions not allowed!!!"
	  << "\n which are: " <<  _wignerDRefKey << " and: " << ref 
	  << "\n clone particle(s) for at least one decay chain" << endmsg;
    exit(1);
  }
  _wignerDRefKey=ref;

  if (!_daughter1IsStable){
    _absDecDaughter1->setWigDRefKey(_refKey);
  }
  if (!_daughter2IsStable){
    _absDecDaughter2->setWigDRefKey(_refKey);
  }
  _nameId=IdStringMapRegistry::instance()->stringId(_name);
  _wignerDId=IdStringMapRegistry::instance()->stringId(_wignerDKey);
  std::string wignerDqNormKey=_wignerDKey+"qNorm";
  std::string keyForMassList="Mass";
  _wignerDqId=IdStringMapRegistry::instance()->keyStringId(keyForMassList,_wignerDKey); 
  _wignerDqNormId = IdStringMapRegistry::instance()->keyStringId(keyForMassList, wignerDqNormKey);
  _wigDWigDRefId=IdStringMapRegistry::instance()->stringStringId(_wignerDKey, _wignerDRefKey);
  NoticeMsg << "wigDWigDRefId = " << _wigDWigDRefId << endmsg;

  std::pair<std::string, std::string > wigDWigDRefPair =
    IdStringMapRegistry::instance()->stringPair(_wigDWigDRefId);
  InfoMsg << "name of amplitude: " << _name << "\tnameId: " << _nameId
	  <<"\n_wigDWigDRefId = " << _wigDWigDRefId << "\t_wignerDId = " << _wignerDId 
	  << "\nwith wignerDKey: " << wigDWigDRefPair.first << "\twignerDRefKey: " 
	  << wigDWigDRefPair.second << endmsg;   
}

 void AbsDecay::setMassParKey(const std::string& newKey){
   _massParamKey = newKey;
   _massParamId = IdStringMapRegistry::instance()->keyStringId("grandMaAndMassParKey", _massParamKey);
}

