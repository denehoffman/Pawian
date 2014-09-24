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
#include "PwaUtils/FlatteDynamics.hh"
#include "PwaUtils/KPiSWaveIso12Dynamics.hh"
#include "PwaUtils/KPiSWaveIso32Dynamics.hh"
#include "PwaUtils/PiPiSWaveASDynamics.hh"
#include "PwaUtils/KMatrixDynamics.hh"
#include "PwaUtils/VoigtDynamics.hh"
#include "PwaUtils/K0star1430LassDynamics.hh"
#include "PwaUtils/BlattWBarrierDynamics.hh"

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
  if(0==theDec->motherPart()) theName=theDec->motherIGJPC()->jpcname();
  //  std::string theName=theDec->motherIGJPC()->jpcname();

  std::string dynType=theDec->dynType();

  std::shared_ptr<AbsDynamics> result;

  std::map<std::string, std::shared_ptr<AbsDynamics> >::iterator it = _dynMap.find(theName);
  if (it !=_dynMap.end()){
    result=it->second;
    Info << "Dynamics\t" << theName << "\tfound!!!" << endmsg;
  }
  else{
    std::vector<Particle*> fsParticles=theDec->finalStateParticles();

    if(theDec->hasMother()){

      if(theDec->dynType()=="BreitWigner")
    	result= std::shared_ptr<AbsDynamics>(new BreitWignerDynamics(theName, fsParticles, theDec->motherPart()));
      else if(theDec->dynType()=="BreitWignerRel")
    	result= std::shared_ptr<AbsDynamics>(new BreitWignerRelDynamics(theName, fsParticles, theDec->motherPart(), theDec->massSumFsParticlesDec1(), theDec->massSumFsParticlesDec2() ));
      else if(theDec->dynType()=="BreitWignerBlattWRel")
	result= std::shared_ptr<AbsDynamics>(new BreitWignerBlattWRelDynamics(theName, fsParticles, theDec->motherPart(), theDec->massSumFsParticlesDec1(), theDec->massSumFsParticlesDec2(), theDec->barrierqR()));
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
      else if(theDec->dynType()=="BlattWBarrier") 
	result= std::shared_ptr<AbsDynamics>(new BlattWBarrierDynamics(theName, fsParticles, theDec->motherPart(), theDec->massSumFsParticlesDec1(), theDec->massSumFsParticlesDec2(), theDec->wignerDKey(), theDec->barrierqR()));


      else if(theDec->dynType()=="WoDynamics") result= std::shared_ptr<AbsDynamics>(new WoDynamics(theName, fsParticles, theDec->motherPart()));
      else{
    	Alert << "Dyn type:\t" << theDec->dynType() << "\tdoes not exist" << endmsg;
    	exit(1);
      }
    }
    else{ // has no mother
      if(theDec->dynType()=="WoDynamics") result= std::shared_ptr<AbsDynamics>(new WoDynamics(theName, fsParticles, theDec->motherPart()));
      else{
	Alert << "no mother resonance; can not add dynamis"
	      << "\nDyn type:\t" << theDec->dynType()
	      << endmsg;
	exit(1);
      }
    }
    result->setMassKey(theDec->massParKey());
    _dynMap[theName]=result;
    _dynVec.push_back(result);
  }

  return result;
}


