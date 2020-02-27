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
// Copyright 2017 Bertram Kopf, 2020 Malte Albrecht


#include "ConfigParser/pipiScatteringParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


pipiScatteringParser::pipiScatteringParser(int argc,char **argv):
  ParserBase(argc,argv)
  ,_pathToKMatrixCompareFile("")
  ,_baseNameFVector("")
  ,_baseNameFVectorCompare("")
  ,_minImagMass(-0.4)
  ,_maxImagMass(-0.01)
  ,_minRealMass(0.5)
  ,_maxRealMass(2.0)
  ,_startRealMass(1.2)
  ,_startImagMass(0.15)
  ,_sheet("nn")
{
  po::options_description common("Common Options");
  common.add_options()
    ("pathToKMatrixCompareFile", po::value<std::string>(&_pathToKMatrixCompareFile)->default_value(_pathToKMatrixCompareFile),"path to another KMatrix-cfg file for comparison purposes")
    ("baseNameFVector",po::value<std::string>(&_baseNameFVector)->default_value(_baseNameFVector),"base name of F-vector for comparison")
    ("baseNameFVectorCompare",po::value<std::string>(&_baseNameFVectorCompare)->default_value(_baseNameFVectorCompare),"base name of second F-vector for comparison")
  ("minImagMass",po::value<double>(&_minImagMass)->default_value(_minImagMass),"minimal imaginary mass for Rieman scan; must be nagative")
  ("maxImagMass",po::value<double>(&_maxImagMass)->default_value(_maxImagMass),"maximal imaginary mass for Rieman scan; must be nagative")
  ("minRealMass",po::value<double>(&_minRealMass)->default_value(_minRealMass),"minimal real mass for Rieman scan")
  ("maxRealMass",po::value<double>(&_maxRealMass)->default_value(_maxRealMass),"maximal real mass for Rieman scan")
   ("startRealMass",po::value<double>(&_startRealMass)->default_value(_startRealMass),"startimg value of real mass for Rieman scan")
   ("startImagMass",po::value<double>(&_startImagMass)->default_value(_startImagMass),"starting value of imaginary mass for Rieman scan")
    ("sheet", po::value<std::string>(&_sheet)->default_value(_sheet),"specification of the Rieman sheet for the scan")
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
  if (_minImagMass>0) _minImagMass = -_minImagMass;
  if (_maxImagMass>0) _maxImagMass = -_maxImagMass;
  
  std::cout << "pathToKMatrixCompareFile = " << _pathToKMatrixCompareFile << std::endl;
  std::cout << "baseNameFVector = " << _baseNameFVector << std::endl;
  std::cout << "baseNameFVectorCompare = " << _baseNameFVectorCompare << std::endl;
  std::cout << "minImagMass = " << _minImagMass << std::endl;
  std::cout << "maxImagMass = " << _maxImagMass << std::endl;
  std::cout << "minRealMass = " << _minRealMass << std::endl;
  std::cout << "maxRealMass = " << _maxRealMass << std::endl;
  std::cout << "startRealMass = " << _startRealMass << std::endl;
  std::cout << "startImagMass = " << _startImagMass << std::endl;
  std::cout << "sheet = " << _sheet << std::endl;
  std::cout << std::endl;

  return true;
}
