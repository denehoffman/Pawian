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

// IdStringMapRegistry class definition file. -*- C++ -*-
// Copyright 2016 Bertram Kopf

#pragma once

#include <iostream>
#include <utility>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <memory>

class IdStringMapRegistry{

public:

  /** Destructor */
  virtual ~IdStringMapRegistry();

  static IdStringMapRegistry* instance();

  unsigned short stringStringId(std::string&, std::string&);
  std::pair<std::string, std::string > stringPair(unsigned short id);

  unsigned short stringId(std::string&);
  std::string getString(unsigned short id);

protected:
 ///Constructor
  IdStringMapRegistry();
  static IdStringMapRegistry* _instance;

private:

  std::map<unsigned short, std::pair<std::string, std::string > > _stringPairMap;
  unsigned short _stringPairMapCounter;

  std::map<unsigned short, std::string > _stringMap;
  unsigned short _stringMapCounter;
};



