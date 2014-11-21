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

// epemParser class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf


#include "ConfigParser/epemParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


epemParser::epemParser(int argc,char **argv):
  ParserBase(argc,argv)
  ,_cmsMass(3.096916)
  ,_productionTensorRadType("Default")
{
  po::options_description common("Common Options");
  common.add_options()
     ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ("spinDensity", po::value< vector<string> >(&_spinDensity), "particles for spin density matrix calculation") 
    ("cmsMass", po::value<double>(&_cmsMass), "CMS mass")
    ("productionTensorRadType", po::value< string >(&_productionTensorRadType), "chosen tensor formalism for describing radiative decays in the production: Zou, Koch")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool epemParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);
  
  std::cout << "\ncms mass:\t" << _cmsMass << std::endl;
  
  std::vector<std::string>::const_iterator it;
  
  std::cout << "\nspin density matrix calculation for particles" << std::endl;
  for (it=_spinDensity.begin(); it!=_spinDensity.end(); ++it){
    std::cout << (*it) << "\n";
  }
  
  std::cout << "\nproductionTensorRadType: " << _productionTensorRadType << std::endl;  
  std::cout << std::endl;
  
  return true;
}
