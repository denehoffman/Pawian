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

// IsobarHeliDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/IsobarHeliDecay.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

IsobarHeliDecay::IsobarHeliDecay(Particle* mother, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv) :
  AbsDecay(mother, daughter1, daughter2, theEnv)
{
}

IsobarHeliDecay::IsobarHeliDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, AbsEnv* theEnv, std::string motherName) :
  AbsDecay(motherIGJPCPtr, daughter1, daughter2, theEnv, motherName)
{
}

void IsobarHeliDecay::extractStates(){
  if (_useIsospin){
    Spin currentGParity=_motherIGJPCPtr->G;
    if(currentGParity==-1 || currentGParity==1) validJPClamlam( _motherIGJPCPtr, _daughter1, _daughter2, _JPClamlamDecAmps, true, _gParity, true);
  }
  else validJPClamlam( _motherIGJPCPtr, _daughter1, _daughter2, _JPClamlamDecAmps); 
}

IsobarHeliDecay::~IsobarHeliDecay(){
}

void IsobarHeliDecay::print(std::ostream& os) const{
  os << "\nJPClamlam amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";
  
  std::vector< std::shared_ptr<const JPClamlam> >::const_iterator it;
  for (it = _JPClamlamDecAmps.begin(); it!= _JPClamlamDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }

  AbsDecay::print(os);  
  os << "\n";
}
