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

// Bertram Kopf (RUB)

#include "ConfigParser/KMatrixParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


KMatrixParser::KMatrixParser(std::string& path)
      : _configFile(path)
      , _noOfChannels(0)
      , _noOfPoles(0)
      ,_config(new po::options_description("Configuration file options"))

     {
    _config->add_options()
      ("configFile",po::value<std::string>(&_configFile),
     	    "The name of the configuration file holding further configuration options")
      ("noOfChannels",po::value<unsigned int>(&_noOfChannels),  "number of channels")
      ("noOfPoles",po::value<unsigned int>(&_noOfPoles),  "number of poles")
      ("gFactor",po::value< vector<string> >(&_gFactors),  "g-factors for the individual channels")
      ("pole",po::value< vector<string> >(&_poles),  "pole: name and mass")
      ("projection",po::value< string>(&_projection),  "projection of the P-vector via pair of decay particles")
       ;


    parseCommandLine();
  }

/************************************************************************************************/
/************************************************************************************************/
/**
  * A function that parses the command line for all required parameters
  */

bool KMatrixParser::parseCommandLine()
{
    po::options_description config_file_options;
    config_file_options.add(*_config);

    std::ifstream ifs(_configFile.c_str());
    if(!ifs.good())
      {
    	stringstream strError;
    	strError << "Error accessing configuratiocommonn file " << _configFile;
    	std::cout << config_file_options << endl;
	exit(0);
      }

    po::variables_map vm;
    store(po::parse_config_file(ifs, config_file_options), vm);
    po::notify(vm);

    std::cout << "\nRunning with the following options using " << _configFile << ":\n\n"
	      << "number of channels: " << _noOfChannels  << "\n\n"
	      << "number of poles: " << _noOfPoles  << "\n\n"
	      << endl;
    
    std::cout << "g-factors are defined as follows:" << std::endl;
    std::vector<std::string>::const_iterator it;
    for(it=_gFactors.begin(); it!=_gFactors.end(); ++it){
      std::cout << (*it) << "\n";
    }

    std::cout << "\npoles: name and mass" << std::endl;
    for(it=_poles.begin(); it!=_poles.end(); ++it){
      std::cout << (*it) << "\n";
    }

    std::cout << "\nprojection: " << _projection << std::endl;
 
  return true;
}

