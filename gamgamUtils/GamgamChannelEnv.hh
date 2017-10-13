//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// GamgamChannelEnv class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include "PwaUtils/AbsChannelEnv.hh"


class gamgamParser;
class gamgamReaction;


class GamgamChannelEnv : public AbsChannelEnv
{
public:
   virtual void setupChannel(ChannelID id);
   GamgamChannelEnv(gamgamParser* theParser);

   std::shared_ptr<gamgamReaction> reaction() {return _gamgamReaction;}
   const double cmsMass() {return _cmsMass;}
   virtual const std::string  channelTypeName() {return "gamgam";}

protected:

private:
   gamgamParser* _theGamGamParser;
   double _cmsMass;
  unsigned _jmax;
   std::shared_ptr<gamgamReaction> _gamgamReaction;
};
