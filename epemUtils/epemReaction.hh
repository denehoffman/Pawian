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

// epemReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <memory>

#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class IsobarLSDecay;
class IsobarHeliDecay;

class epemReaction {

public:
  epemReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, ChannelID channelID);

  virtual ~epemReaction();

  virtual void print(std::ostream& os) const;
  std::vector< std::shared_ptr<IsobarLSDecay> >& productionCanoDecays() {return _prodCanoDecs;}
  std::vector< std::shared_ptr<IsobarHeliDecay> >& productionHeliDecays() {return _prodHeliDecs;}

protected:

private:
  ChannelID _channelID;
  std::shared_ptr<const IGJPC> _epemIGJPC;

  std::vector< std::shared_ptr<IsobarLSDecay> > _prodCanoDecs;
  std::vector< std::shared_ptr<IsobarHeliDecay> > _prodHeliDecs;

};
