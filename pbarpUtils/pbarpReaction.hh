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

// pbarpReaction class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class Particle;
class IsobarLSDecay;
class IsobarHeliDecay;
class IsobarTensorDecay;

class pbarpReaction {

public:
  pbarpReaction(std::vector<std::pair<Particle*, Particle*> >& prodPairs, int lmax=3);

  virtual ~pbarpReaction();

  virtual void print(std::ostream& os) const;
  std::vector< boost::shared_ptr<IsobarLSDecay> >& productionDecays() {return _prodDecs;}
  std::vector< boost::shared_ptr<IsobarHeliDecay> >& productionHeliDecays() {return _prodHeliDecs;}
  std::vector< boost::shared_ptr<IsobarTensorDecay> >& productionTensorDecays() {return _prodTensorDecs;}  
  std::vector< boost::shared_ptr<const jpcRes> >& jpcStates() {return _pbarpJPCs;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsStates() {return _pbarpJPCLSs;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsSingletStates() {return _pbarpJPCLSsinglet;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsTriplet0States() {return _pbarpJPCLStriplet0;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsTripletp1States() {return _pbarpJPCLStripletp1;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsTripletm1States() {return _pbarpJPCLStripletm1;}

protected:

private:
  int _lmax;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpJPCs;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLSs;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLSsinglet;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLStriplet0;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLStripletp1;
  std::vector< boost::shared_ptr<const JPCLS> > _pbarpJPCLStripletm1;

  std::vector< boost::shared_ptr<IsobarLSDecay> > _prodDecs;
  std::vector< boost::shared_ptr<IsobarHeliDecay> > _prodHeliDecs;
  std::vector< boost::shared_ptr<IsobarTensorDecay> > _prodTensorDecs;

  std::map< boost::shared_ptr<const jpcRes>, std::vector< boost::shared_ptr<const JPCLS> >, pawian::Collection::SharedPtrLess> _jpcToJPCLSMap;

  bool CheckJPCLSForParticle(std::string& particleName, boost::shared_ptr<const JPCLS> theJPCLS);
};
