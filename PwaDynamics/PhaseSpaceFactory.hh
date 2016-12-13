//************************************************************************//
//									  //
//  Copyright 2016 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// PhaseSpaceFactory class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#pragma once

#include <iostream>
#include <string>

// #include <cassert>
#include <memory>
#include "PwaDynamics/AbsPhaseSpace.hh"

class AbsPhaseSpace;

class PhaseSpaceFactory{

public:

  /** Destructor */
  virtual ~PhaseSpaceFactory();

  static PhaseSpaceFactory* instance();
  std::shared_ptr<AbsPhaseSpace> getPhpPointer(std::string type, std::vector<double> masses);

protected:
 ///Constructor 
  PhaseSpaceFactory();
  static PhaseSpaceFactory* _instance;

private:
};



