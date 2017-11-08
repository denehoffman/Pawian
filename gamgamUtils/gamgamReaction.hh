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

// gamgamReaction class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

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
#include "PwaUtils/FormationDecay.hh"
#include "Utils/PawianCollectionUtils.hh"


class Particle;
class IsobarHeliDecay;
class ProdChannelInfo;
class gamgamStates;

class gamgamReaction {

public:
  gamgamReaction(std::vector<std::shared_ptr<ProdChannelInfo> > prodChannelInfoList, ChannelID channelID, int jmax);

  virtual ~gamgamReaction();

  virtual void print(std::ostream& os) const;
  std::shared_ptr<gamgamStates> GamGamStates() {return _gamgamStates;}
  std::vector< std::shared_ptr<FormationDecay> >& formationDecays() {return _prodFormationDecs;}
  std::vector< std::shared_ptr<IsobarHeliDecay> >& productionHeliDecays() {return _prodHeliDecs;}
  std::shared_ptr<FormationDecay> motherProdDec(){return _motherProdDec;}

protected:

private:
  ChannelID _channelID;
  unsigned int _jmax;
  std::shared_ptr<gamgamStates> _gamgamStates;
  std::vector< std::shared_ptr<FormationDecay> > _prodFormationDecs;
  std::vector< std::shared_ptr<IsobarHeliDecay> > _prodHeliDecs;
  std::shared_ptr<FormationDecay> _motherProdDec;
};
