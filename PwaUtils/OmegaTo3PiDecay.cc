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

// OmegaTo3PiDecay class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <algorithm>

#include "PwaUtils/OmegaTo3PiDecay.hh"
#include "PwaUtils/AbsEnv.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "Utils/FunctionUtils.hh"
#include "PwaUtils/KinUtils.hh"
#include "PwaUtils/EvtDataBaseList.hh"

OmegaTo3PiDecay::OmegaTo3PiDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, AbsEnv* theEnv) :
  AbsDecay(mother, daughter1, daughter2, theEnv)
  ,_daughter3(daughter3)
{
  _finalStateParticles.push_back(daughter3);
  _finalStateParticlesDaughter3.push_back(daughter3);

  pawian::Collection::PtrLess thePtrLess;
  std::sort(_finalStateParticles.begin(), _finalStateParticles.end(), thePtrLess);

  _name+="_"+daughter3->name();
  _fitParamSuffix=_name;
  
  //check correct quantum numbers
  //...
  Spin validL=1;
  Spin validS=0;
  std::shared_ptr<const JPCLS> theValidJPCLS(new JPCLS(_motherJPCPtr, validL, validS));
  _JPCLSDecAmps.push_back(theValidJPCLS);

}

OmegaTo3PiDecay::~OmegaTo3PiDecay(){
}

void OmegaTo3PiDecay::fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtData* evtData){
  return;
}

void OmegaTo3PiDecay::print(std::ostream& os) const{
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
