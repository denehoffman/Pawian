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

// EpemChannelEnv class definition file. -*- C++ -*-
// Copyright 2013 Julian Pychy

#pragma once

#include "PwaUtils/AbsChannelEnv.hh"


class epemParser;
class epemReaction;


class EpemChannelEnv : public AbsChannelEnv
{
public:
   void setup(ChannelID id);
   EpemChannelEnv(epemParser* theParser);

   std::shared_ptr<epemReaction> reaction() {return _epemReaction;}
   std::vector<std::string>& spinDensityNames(){ return _spinDensity;}
   const double cmsMass() {return _cmsMass;}
  virtual const std::string  channelTypeName() {return "epem";}
   void CreateHistInstance(std::shared_ptr<AbsLh> theLh, fitParams& theFitParams, std::string additionalSuffix);

private:
   epemParser* _theParser;
   double _cmsMass;
   std::shared_ptr<epemReaction> _epemReaction;
   std::vector<std::string> _spinDensity;
};
