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

// OmegaTo3PiLSDecay class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/DataUtils.hh"
#include "Utils/PawianCollectionUtils.hh"

class Particle;
class EvtData;
class AbsEnv;

class OmegaTo3PiLSDecay : public AbsDecay{

public:
  OmegaTo3PiLSDecay(Particle* mother, Particle* daughter1, Particle* daughter2, Particle* daughter3, AbsEnv* theEnv);
  virtual ~OmegaTo3PiLSDecay();
  //  virtual OmegaTo3PiLSDecay* clone_() const = 0;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSAmps(){ return _JPCLSDecAmps;}
  virtual void fillWignerDs(std::map<std::string , Vector4<double> >& fsMap, EvtData* evtData);
  virtual const std::string lambdaDecKey() {return _lambdaDecKey;}
  virtual void print(std::ostream& os) const;
  Particle* daughter3Part() {return _daughter3;}
  virtual std::string type() {return "OmegaTo3PiLSDecay";} 

protected:
  std::string _lambdaDecKey;
  Particle* _daughter3;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSDecAmps;
  std::vector<Particle*> _finalStateParticlesDaughter3;

};
