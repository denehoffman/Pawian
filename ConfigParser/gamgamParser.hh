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

// gamgamParser class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#pragma once

#include "ConfigParser/ParserBase.hh"
// Boost headers go here



class gamgamParser : public ParserBase 
{

  public:

  gamgamParser(int argc,char **argv);
  virtual ~gamgamParser(){;}

  const double cmsMass() {return _cmsMass;}
  const unsigned int jMax() {return _jMax;}

protected:
  virtual bool parseCommandLine(int argc,char **argv); 

  double _cmsMass;
  unsigned int _jMax;
};



