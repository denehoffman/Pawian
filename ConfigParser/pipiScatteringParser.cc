//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// pipiScatteringParser class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf


#include "ConfigParser/pipiScatteringParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


pipiScatteringParser::pipiScatteringParser(int argc,char **argv):
  ParserBase(argc,argv)
{
  po::options_description common("Common Options");
  common.add_options()
    ("scatterDummy", po::value<float>(&_scatterDummy)->default_value(_scatterDummy),"antiproton momentum")
    ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool pipiScatteringParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);

  std::cout << "scatterDummy = " << _scatterDummy << std::endl;

  std::cout << std::endl;

  return true;
}
