//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// IsobarTensorReducedRadDecay class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/AbsDynamics.hh"
#include "PwaUtils/IsobarTensorReducedRadDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/DynRegistry.hh"
#include "ConfigParser/ParserBase.hh"

IsobarTensorReducedRadDecay::IsobarTensorReducedRadDecay(Particle* mother, Particle* daughter1_gamma, Particle* daughter2, ChannelID channelID) :
  IsobarTensorDecay(mother, daughter1_gamma, daughter2, channelID)
  ,_XisEven(false)
  ,_noOfAmps(0)
{
  if(mother->twoIso() !=0 || mother->theParity() != -1 || mother->twoJ() != 2 || mother->theParity() != -1 ||  mother->theCParity() !=-1 || mother->mass() < 0.01){
    Alert << "so far only mother particle with quantum number combinations JPC=1-- or its masses abov 0.01 GeV/c2 supported!!!" << endmsg;
    exit(0); 
  }
 if(daughter1_gamma->twoJ() != 2 || daughter1_gamma->theParity() != -1 ||  daughter1_gamma->theCParity() || daughter1_gamma->mass() < 1.e-6){
   Alert << "daughter particle 1 is not a photon" << endmsg;
   exit(0); 
 }

 if(daughter2->twoJ()%4 == 0) _XisEven=true;
 if(_XisEven){
   Alert << "decay particle with even J not supported so far!!!" << endmsg;
   exit(0);
 }

 if(daughter2->twoJ() != 2 || daughter2->theParity() != 1){
   Alert << "decay particle with even JP!=1+ not supported so far!!!" << endmsg;
   exit(0);
 }
 
 if(daughter2->twoJ() == 2 && daughter2->theParity() == 1){
     _noOfAmps=2;
   }
}

IsobarTensorReducedRadDecay::IsobarTensorReducedRadDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1_gamma, Particle* daughter2, ChannelID channelID, std::string motherName) :
  IsobarTensorDecay(motherIGJPCPtr, daughter1_gamma, daughter2, channelID, motherName)
  ,_XisEven(false)
  ,_noOfAmps(0)
{
  if(motherIGJPCPtr->I !=0 || motherIGJPCPtr->G != -1 || motherIGJPCPtr->J != 1 || motherIGJPCPtr->P != -1 ||  motherIGJPCPtr->C !=-1){
    Alert << "so far only mother particle with quantum number combinations JPC=1-- !!!" << endmsg;
    exit(0); 
  }
  if(daughter1_gamma->twoJ() != 2 || daughter1_gamma->theParity() != -1 ||  daughter1_gamma->theCParity() || daughter1_gamma->mass() < 1.e-6){
    Alert << "daughter particle 1 is not a photon" << endmsg;
    exit(0); 
  }
  
  if(daughter2->twoJ()%4 == 0) _XisEven=true;
  if(_XisEven){
    Alert << "decay particle with even J not supported so far!!!" << endmsg;
    exit(0);
  }
  
  if(daughter2->twoJ() != 2 || daughter2->theParity() != 1){
    Alert << "decay particle with even JP!=1+ not supported so far!!!" << endmsg;
    exit(0);
  }
  
  if(daughter2->twoJ() == 2 && daughter2->theParity() == 1){
    _noOfAmps=2;
  }
}

IsobarTensorReducedRadDecay::~IsobarTensorReducedRadDecay(){
}

void IsobarTensorReducedRadDecay::print(std::ostream& os) const{
  //dummy
}

void IsobarTensorReducedRadDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData, std::string& refKey){
  //dummy
}

