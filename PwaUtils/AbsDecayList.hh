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

// AbsDecayList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <memory>

#include "PwaUtils/DataUtils.hh"

class Particle;
class AbsDecay;

class AbsDecayList {

public:

  AbsDecayList();
  ~AbsDecayList();

  void addDecay(std::shared_ptr<AbsDecay> theIsoDec);
  std::shared_ptr<AbsDecay> decay(Particle* mother);
  std::shared_ptr<AbsDecay> decay(const std::string& name);
  void replaceSuffix(const std::string& oldPart, const std::string& newPart);
  void replaceMassKey(const std::string& oldPart, const std::string& newPart);
  std::vector<std::shared_ptr<AbsDecay> >& getList() {return _absDecList;}    
protected:  

  std::vector<std::shared_ptr<AbsDecay> > _absDecList;  
};
