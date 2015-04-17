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

// IsobarLSDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/IsobarLSDecay.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

IsobarLSDecay::IsobarLSDecay(Particle* mother, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string typeName) :
  AbsDecay(mother, daughter1, daughter2, channelID, typeName)
{
}

IsobarLSDecay::IsobarLSDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, Particle* daughter2, ChannelID channelID, std::string motherName, std::string typeName) :
  AbsDecay(motherIGJPCPtr, daughter1, daughter2, motherName, channelID, typeName)
{
}

void IsobarLSDecay::extractStates(){
  extractLmin();
  if (_useIsospin){
    Spin currentGParity=_motherIGJPCPtr->G;
    int daughter1GParity=_daughter1->theGParity();
    int daughter2GParity=_daughter2->theGParity();

    if( fabs(currentGParity)==1 && fabs(daughter1GParity)==1 && fabs(daughter2GParity)==1){
      validJPCLS( _motherIGJPCPtr, _daughter1, _daughter2, _JPCLSDecAmps, true, _gParity, true);
      validLS( _motherIGJPCPtr, _daughter1, _daughter2, _LSDecAmps, true, _gParity, true);
    }
    else{
      validJPCLS( _motherIGJPCPtr, _daughter1, _daughter2, _JPCLSDecAmps);
      validLS( _motherIGJPCPtr, _daughter1, _daughter2, _LSDecAmps);
    }
  }
  else{
    validJPCLS( _motherIGJPCPtr, _daughter1, _daughter2, _JPCLSDecAmps);
    validLS( _motherIGJPCPtr, _daughter1, _daughter2, _LSDecAmps);
  }

  if( 0==_JPCLSDecAmps.size()) Info << "_JPCLSDecAmps.size()==0 for " << name() << endmsg;
  print(std::cout);
}

IsobarLSDecay::~IsobarLSDecay(){
}

void IsobarLSDecay::print(std::ostream& os) const{
  os << "\nJPCLS amplitudes for decay\t" << _name << ":\n";
  os << "suffix for fit parameter name:\t" << _fitParamSuffix << "\n";

  std::vector< std::shared_ptr<const JPCLS> >::const_iterator it;
  for (it = _JPCLSDecAmps.begin(); it!= _JPCLSDecAmps.end(); ++it){
    (*it)->print(os);
    os << "\n";
  }

  AbsDecay::print(os);
  os << "\n";
}
