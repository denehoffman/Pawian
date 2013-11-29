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

// gammapReaction class definition file. -*- C++ -*-
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
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class Particle;
class IsobarLSDecay;
class IsobarHeliDecay;
class IsobarTensorDecay;

class gammapReaction {

public:
  gammapReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, ChannelID channelID, int lmax=3);

  virtual ~gammapReaction();

  virtual void print(std::ostream& os) const;
  std::vector< std::shared_ptr<IsobarLSDecay> >& productionDecays() {return _prodLSDecs;}
  std::vector< std::shared_ptr<IsobarHeliDecay> >& productionHeliDecays() {return _prodHeliDecs;}
  std::vector< std::shared_ptr<IsobarTensorDecay> >& productionTensorDecays() {return _prodTensorDecs;}
  std::vector< std::shared_ptr<const IGJPC> >& igjpcStates() {return _gammapIGJPCs;}
  std::vector< std::shared_ptr<const JPCLS> >& JPCljStates() {return _gammapJPCljs;}
  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const JPCLS> >, pawian::Collection::SharedPtrLess>& jpcToJPCljMap(){return _jpcToJPCljMap;}
  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const IGJPC> >, pawian::Collection::SharedPtrLess>& jpcToIGJPCMap(){return _jpcToIGJPCMap;}
  
protected:

private:
  ChannelID _channelID;
  int _lmax;
  std::vector< std::shared_ptr<const IGJPC> > _gammapIGJPCs;
  std::vector< std::shared_ptr<const JPCLS> > _gammapJPCljs;

  std::vector< std::shared_ptr<IsobarLSDecay> > _prodLSDecs;
  std::vector< std::shared_ptr<IsobarHeliDecay> > _prodHeliDecs;
  std::vector< std::shared_ptr<IsobarTensorDecay> > _prodTensorDecs;

  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const JPCLS> >, pawian::Collection::SharedPtrLess> _jpcToJPCljMap;
  std::map< std::shared_ptr<const jpcRes>, std::vector< std::shared_ptr<const IGJPC> >, pawian::Collection::SharedPtrLess> _jpcToIGJPCMap;

  bool CheckJPCLSForParticle(std::string& particleName, std::shared_ptr<const JPCLS> theJPCLS);
};
