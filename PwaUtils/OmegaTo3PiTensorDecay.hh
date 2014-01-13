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

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <memory>

#include "PwaUtils/OmegaTo3PiDecay.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;

class OmegaTo3PiTensorDecay : public OmegaTo3PiDecay{

public:
  OmegaTo3PiTensorDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, ChannelID channelID);
  virtual ~OmegaTo3PiTensorDecay();

  virtual std::string type() {return "OmegaTo3PiTensorDecay";}

  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData);

protected:
};
