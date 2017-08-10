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


#include "ConfigParser/gamgamParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


gamgamParser::gamgamParser(int argc,char **argv):
  ParserBase(argc,argv)
  ,_cmsMass(3.096916)
  ,_jMax(2)
{
  po::options_description common("Common Options");
  common.add_options()
     ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ("cmsMass", po::value<double>(&_cmsMass), "CMS mass")
    ("jMax", po::value<unsigned int>(&_jMax), "maximal spin of the formed resonances")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool gamgamParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);
  
  std::cout << "\ncms mass:\t" << _cmsMass << std::endl;
  std::cout << "\njMax:\t" << _jMax << std::endl;
  
  std::cout << std::endl;
  
  return true;
}
