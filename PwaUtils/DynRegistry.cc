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

// DynRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/DynRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/BreitWignerDynamics.hh"
#include "PwaUtils/BreitWignerRelDynamics.hh"
#include "PwaUtils/BreitWignerBlattWRelDynamics.hh"
#include "PwaUtils/BreitWignerBlattWTensorRelDynamics.hh"
#include "PwaUtils/FlatteDynamics.hh"
#include "PwaUtils/KPiSWaveIso12Dynamics.hh"
#include "PwaUtils/KPiSWaveIso32Dynamics.hh"
#include "PwaUtils/PiPiSWaveASDynamics.hh"
#include "PwaUtils/KMatrixDynamics.hh"
#include "PwaUtils/VoigtDynamics.hh"
#include "PwaUtils/K0star1430LassDynamics.hh"
#include "PwaUtils/BlattWBarrierDynamics.hh"
#include "PwaUtils/RadM1Dynamics.hh"
#include "PwaUtils/BlattWBarrierTensorDynamics.hh"

#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/WoDynamics.hh"

#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLogger.hh"


DynRegistry* DynRegistry::_instance=0;

DynRegistry* DynRegistry::instance()
{
  if (0==_instance) _instance = new DynRegistry();
  return _instance;
}

DynRegistry::DynRegistry()
{
}

DynRegistry::~DynRegistry()
{
}

std::shared_ptr<AbsDynamics> DynRegistry::getDynamics(std::shared_ptr<AbsDecay> theDec){

  std::string theName=theDec->name();
  //  if(theDec->dynType()=="BlattWBarrier") theName=theDec->wignerDKey();
  if(0==theDec->motherPart()){
    if(theDec->dynType()=="WoDynamics") theName=theDec->motherIGJPC()->jpcname();
    else if (theDec->dynType()!="BlattWBarrier" && theDec->dynType()!="BlattWBarrierTensor" && theDec->dynType()!="RadM1"){
      Alert << "dynamics with type " << theDec->dynType() << " is not allowed for annihilation amplitude!!!" << endmsg;
      exit(0); 
    }
  }

  std::string dynType=theDec->dynType();

  std::shared_ptr<AbsDynamics> result;

  ChannelID currentChannelId = theDec->channelId();
  std::map<std::string, std::shared_ptr<AbsDynamics> > currentDynMap;

  std::map<ChannelID, std::map<std::string, std::shared_ptr<AbsDynamics> > >::iterator itChannelId =_dynMapChannel.find(currentChannelId);

  if (itChannelId !=_dynMapChannel.end()){
    currentDynMap=itChannelId->second;
  }
  else _dynMapChannel[currentChannelId]=currentDynMap;

  std::map<std::string, std::shared_ptr<AbsDynamics> >::iterator it = currentDynMap.find(theName);
  if (it !=currentDynMap.end()){
    result=it->second;
    Info << "Dynamics\t" << theName << "\tfound!!!" << endmsg;
  }
  else{
    std::vector<Particle*> fsParticles=theDec->finalStateParticles();
    std::vector<Particle*> fsParticlesDaughter1=theDec->finalStateParticlesDaughter1();
    std::vector<Particle*> fsParticlesDaughter2=theDec->finalStateParticlesDaughter2();

    if(theDec->dynType()=="BreitWigner")
      result= std::shared_ptr<AbsDynamics>(new BreitWignerDynamics(theName, fsParticles, theDec->motherPart()));
    else if(theDec->dynType()=="BreitWignerRel")
      result= std::shared_ptr<AbsDynamics>(new BreitWignerRelDynamics(theName, fsParticles, theDec->motherPart(), fsParticlesDaughter1, fsParticlesDaughter2));
    else if(theDec->dynType()=="BreitWignerBlattWRel"){
      if(theDec->isTensorAmp()){
	Alert << "dynamics BreitWignerBlattWRel is not allowed for non tensor amplitudes (amp name: " << theDec->name() << endmsg;
	exit(0);
      }
      result= std::shared_ptr<AbsDynamics>(new BreitWignerBlattWRelDynamics(theName, fsParticles, theDec->motherPart(), fsParticlesDaughter1, fsParticlesDaughter2, theDec->barrierqR()));
    }
    else if(theDec->dynType()=="BreitWignerBlattWTensorRel"){
      if(!theDec->isTensorAmp()){
	Alert << "dynamics BreitWignerBlattWTensorRel is not allowed for tensor amplitudes (amp name: " << theDec->name() << endmsg;
	exit(0);
      }
      result= std::shared_ptr<AbsDynamics>(new BreitWignerBlattWTensorRelDynamics(theName, fsParticles, theDec->motherPart(), fsParticlesDaughter1, fsParticlesDaughter2, theDec->barrierqR()));
    }
    else if(theDec->dynType()=="KMatrix"){
      std::string pathToConfigFile=theDec->pathToConfigParser();
      result= std::shared_ptr<AbsDynamics>(new KMatrixDynamics(theName, fsParticles, theDec->motherPart(), pathToConfigFile)); 
    }
    else if(theDec->dynType()=="Flatte")
      result= std::shared_ptr<AbsDynamics>(new FlatteDynamics(theName, fsParticles, theDec->motherPart(), theDec->firstDecayChannel(), theDec->secondDecayChannel()));
    else if(theDec->dynType()=="KpiSWaveIso12")
      result= std::shared_ptr<AbsDynamics>(new KPiSWaveIso12Dynamics(theName, fsParticles, theDec->motherPart()));
    else if(theDec->dynType()=="KpiSWaveIso32")
      result= std::shared_ptr<AbsDynamics>(new KPiSWaveIso32Dynamics(theName, fsParticles, theDec->motherPart()));
    else if(theDec->dynType()=="PiPiSWaveAS")
      result= std::shared_ptr<AbsDynamics>(new PiPiSWaveASDynamics(theName, fsParticles, theDec->motherPart(), GlobalEnv::instance()->particleTable()));
    else if(theDec->dynType()=="Voigt") 
      result= std::shared_ptr<AbsDynamics>(new VoigtDynamics(theName, fsParticles, theDec->motherPart()));
    else if(theDec->dynType()=="K0star1430Lass") 
      result= std::shared_ptr<AbsDynamics>(new K0star1430LassDynamics(theName, fsParticles, theDec->motherPart()));
    else if(theDec->dynType()=="BlattWBarrier"){
      if(theDec->isTensorAmp()){
	Alert << "dynamics BlattWBarrier is not allowed for tensor amplitudes (amp name: " << theDec->name() << endmsg;
	exit(0);
      }
      //      std::string blattWBarrierName=theDec->wignerDKey();
      std::string blattWBarrierName=theName;
      result= std::shared_ptr<AbsDynamics>(new BlattWBarrierDynamics(blattWBarrierName, fsParticles, theDec->motherPart(), theDec->wignerDKey(), theDec->barrierqR()));
    }
    else if(theDec->dynType()=="BlattWBarrierTensor"){
      if(!theDec->isTensorAmp()){
	Alert << "dynamics BlattWBarrierTensor is not allowed for tensor amplitudes (amp name: " << theDec->name() << endmsg;
	exit(0);
      }
      // std::string blattWBarrierName=theDec->wignerDKey();
      std::string blattWBarrierName=theName;
      result= std::shared_ptr<AbsDynamics>(new BlattWBarrierTensorDynamics(blattWBarrierName, fsParticles, theDec->motherPart(), theDec->wignerDKey(), theDec->barrierqR()));
    }
    else if(theDec->dynType()=="WoDynamics") result= std::shared_ptr<AbsDynamics>(new WoDynamics(theName, fsParticles, theDec->motherPart()));
    else if(theDec->dynType()=="RadM1") {
      result= std::shared_ptr<AbsDynamics>(new RadM1Dynamics(theName, fsParticles, theDec->motherPart(), fsParticlesDaughter1, fsParticlesDaughter2));
    }
    else{
      Alert << "Dyn type:\t" << theDec->dynType() << "\tdoes not exist" << endmsg;
      exit(1);
    }
    
    result->setMassKey(theDec->massParKey());

    if(theDec->isProductionAmp() &&  (theDec->dynType()=="BlattWBarrier" || theDec->dynType()=="BlattWBarrierTensor" || theDec->dynType()=="RadM1")) result->setMassKey(theDec->prodParKey()); 

    result->fillParamNameList();
    Info << "add dynamics for " <<  theName << endmsg;

    _dynMapChannel[currentChannelId][theName]=result;
    _dynVec.push_back(result);
  }

  return result;
}


