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

// OmegaTo3PiTensorDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/OmegaTo3PiTensorDecay.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

OmegaTo3PiTensorDecay::OmegaTo3PiTensorDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, AbsEnv* theEnv) :
  OmegaTo3PiDecay(mother, daughter1, daughter2, daughter3, theEnv)
{
}

OmegaTo3PiTensorDecay::~OmegaTo3PiTensorDecay(){
}

void OmegaTo3PiTensorDecay::fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtData* evtData){
  int evtNo=evtData->evtNo;
  std::map<int, bool>::const_iterator it = _alreadyFilledMap.find(evtNo);
  if(it!=_alreadyFilledMap.end() &&  it->second) return; //already filled

  std::map<std::string , Vector4<double> >::iterator itMap;

  Vector4<double> all4Vec(0.,0.,0.,0.);
  //fill all4Vec
  for(itMap=fsMap.begin(); itMap!=fsMap.end(); ++itMap){
    all4Vec+=itMap->second;
  }

  //fill daughter1 and daughter2 4Vec
  itMap=fsMap.find(_daughter1->name());
  Vector4<double> daughter1_4Vec=itMap->second;
  // daughter1_4Vec.Boost(all4Vec);

  itMap=fsMap.find(_daughter2->name());
  Vector4<double> daughter2_4Vec=itMap->second;
  // daughter2_4Vec.Boost(all4Vec);
  
  itMap=fsMap.find(_daughter3->name());
  Vector4<double> daughter3_4Vec=itMap->second;
  // daughter3_4Vec.Boost(all4Vec);

  Vector4<double> P_3particle_4Vec=daughter1_4Vec+daughter2_4Vec+daughter3_4Vec;

   Vector4<double> daughter1Hel=helicityVec(all4Vec, P_3particle_4Vec, daughter1_4Vec);
   Vector4<double> daughter2Hel=helicityVec(all4Vec, P_3particle_4Vec, daughter2_4Vec);
   Vector4<double> daughter3Hel=helicityVec(all4Vec, P_3particle_4Vec, daughter3_4Vec);
   Vector4<double> mother_4Vec=helicityVec(all4Vec, P_3particle_4Vec, P_3particle_4Vec);

  LeviCivitaTensor eps; 
  PolVector omega; // spin-1 particle is the default constructor
  //  omega.SetP4(P_3particle_4Vec, P_3particle_4Vec.M());
  omega.SetP4(mother_4Vec, mother_4Vec.M());

  for (Spin mz=-1; mz<=1; ++mz){
    Tensor<complex<double> >  ampTensor;
    //    ampTensor = eps|(daughter1_4Vec%daughter2_4Vec%daughter3_4Vec%omega(mz));
 ampTensor = eps|(daughter1Hel%daughter2Hel%daughter3Hel%omega(mz));
 //    DebugMsg << "mz: " << mz << "\t" << ampTensor << endmsg;    
    evtData->ComplexDoubleString["omegTensor"][_motherJPCPtr->J][mz]=(complex<double>) ampTensor(0);
  }
   _alreadyFilledMap[evtNo]=true;
}

