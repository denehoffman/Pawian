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

// DynRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

class AbsDecay;
class AbsDynamics;
class Particle;
class ParticleTable;

class DynRegistry{

public:

  /** Destructor */
  virtual ~DynRegistry();

  static DynRegistry* instance();
  boost::shared_ptr<AbsDynamics> getDynamics(boost::shared_ptr<AbsDecay> theDec);
  std::vector<boost::shared_ptr<AbsDynamics> > getDynVec(){return _dynVec;}

protected:
 ///Constructor 
  DynRegistry();
  static DynRegistry* _instance;

private:

  std::map<std::string, boost::shared_ptr<AbsDynamics> > _dynMap;
  std::vector<boost::shared_ptr<AbsDynamics> > _dynVec;

};



