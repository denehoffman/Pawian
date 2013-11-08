
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

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include <boost/version.hpp>

#if BOOST_VERSION < 103600
#error "Error: Boost should at least have version 1.36 !"
#endif // BOOST_VERSION

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#pragma once


namespace po = boost::program_options;

class KMatrixParser {
public:
  KMatrixParser(std::string& path);
  virtual ~KMatrixParser(){;}
    
  const std::string& getConfigFile() const { return _configFile;}
  const unsigned int noOfChannels() const {return _noOfChannels;}
  const unsigned int noOfPoles() const {return _noOfPoles;}
  const std::vector<std::string>& gFactors() const { return _gFactors; }
  const std::vector<std::string>& poles() const { return _poles; }
  const std::string projection() const { return _projection; }
  const int useBarrierFactors() const {return _useBarrierFactors;}
  const unsigned int orbitalMom() const {return _orbitalMom;}
  const int orderBg() const {return _orderKMatrixBackground;}
  const bool useAdler() const {return _useAdler0;}
  const double s0Adler() const {return _s0Adler;}
  const double snormAdler() const {return _snormAdler;}
 
protected:
  //  virtual bool parseCommandLine(int argc,char **argv); 
  virtual bool parseCommandLine();
  std::string _configFile;
  unsigned int _noOfChannels;
  unsigned int _noOfPoles;
  std::vector<std::string> _gFactors;
  std::vector<std::string> _poles;
  std::string _projection;
  bool _useBarrierFactors;
  unsigned int _orbitalMom;
  int _orderKMatrixBackground;
  bool _useAdler0;
  double _s0Adler;
  double _snormAdler;
 
  po::options_description* _config;
};
