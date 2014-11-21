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

IsobarTensorReducedRadDecay::IsobarTensorReducedRadDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelID) :
  IsobarTensorDecay(mother, daughter1, daughter2, channelID)
{
}

IsobarTensorReducedRadDecay::IsobarTensorReducedRadDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string motherName) :
  IsobarTensorDecay(motherIGJPCPtr, daughter1, daughter2, channelID, motherName)
{
}

IsobarTensorReducedRadDecay::~IsobarTensorReducedRadDecay(){
}

void IsobarTensorReducedRadDecay::print(std::ostream& os) const{
  //dummy
}

void IsobarTensorReducedRadDecay::fillWignerDs(std::map<std::string, Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData, std::string& refKey){
  //dummy
}

