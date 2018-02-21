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

// FormationDecay class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <memory>

#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/AbsChannelEnv.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;

class FormationDecay : public AbsDecay{

public:
  FormationDecay(std::shared_ptr<const IGJPC> motherIGJPCPtr, Particle* daughter1, ChannelID channelID, std::string motherName="gg", std::string typeName="FormationDecay");
  virtual ~FormationDecay();
  //  virtual FormationDecay* clone_() const = 0;

  virtual void print(std::ostream& os) const;

  virtual void extractStates();
  virtual void setDecayLevel(decLevel theLevel);
  virtual void setDecayLevelTree(decLevel theLevel, std::shared_ptr<AbsDecay> motherDecPtr, std::shared_ptr<AbsDecay> prodDecPtr);
  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, Vector4<double>& prodParticle4Vec, EvtData* evtData);
protected:
};
