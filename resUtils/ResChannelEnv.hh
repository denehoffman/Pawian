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

// ResChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy


#pragma once

//#include <iostream>
//#include <vector>
//#include <map>
//#include <vector>
//#include <string>
//#include <sstream>
//#include <memory>

#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"

class resReaction;
class resParser;
class Particle;

class ResChannelEnv : public AbsChannelEnv{

public:
  void setup(ChannelID id);
  ResChannelEnv(resParser* theResParser);

  std::shared_ptr<resReaction> reaction() {return _resReaction;}
  Particle* motherParticle() {return _motherParticle;}
  virtual const std::string  channelTypeName() {return "res";}
  virtual const bool polarizedMother() const { return _polarizedMother;}
  virtual std::shared_ptr<AbsHist> CreateHistInstance(std::string additionalSuffix);

protected:
  resParser* _theResParser;
  Particle* _motherParticle;
  bool _polarizedMother;
  std::shared_ptr<resReaction> _resReaction;
};
