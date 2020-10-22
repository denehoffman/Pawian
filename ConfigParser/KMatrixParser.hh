
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
  const std::string keyName() const { return _keyName; }   
  const unsigned int noOfChannels() const {return _noOfChannels;}
  const unsigned int noOfPoles() const {return _noOfPoles;}
  const std::vector<std::string>& gFactors() const { return _gFactors; }
  const std::vector<std::string>& phpFactors() const { return _phpDescriptions;}
  const std::vector<std::string>& poles() const { return _poles; }
  const std::string projection() const { return _projection; }
  const int useBarrierFactors() const {return _useBarrierFactors;}
  const bool useTruncatedBarrierFactors() const {return _useTruncatedBarrierFactors;}
  const unsigned int orbitalMom() const {return _orbitalMom;}
  const int orderBg() const {return _orderKMatrixBackground;}
  const int orderPVecBg() const {return _orderPVectorBackground;}
  const bool useAdler() const {return _useAdler0;}
  const double s0Adler() const {return _s0Adler;}
  const double snormAdler() const {return _snormAdler;}
  const std::vector<std::string> gFactorFixPoles() {return _gFactorFixPoles;}
  const std::vector<std::string> gFactorFixReleaseGFactor() {return _gFactorFixReleaseGFactor;}
  const std::vector<std::string> gFactorFixKeepGRatio() {return _gFactorFixKeepGRatio;}
  const std::vector<std::string> gFactorFixMassPol() {return _gFactorFixMassPol;}
  const std::vector<std::string> gFactorFixWidthPol() {return _gFactorFixWidthPol;}
  const std::vector<std::string> gFactorFixSeparateScale() {return _gFactorFixSeparateScale;}
  std::map<std::string, double> gFactorMap() {return _gFactorMap;} 
  std::map<std::pair<std::string, std::string>, std::string> phpDescriptionMap()& {return _phpDescriptionMap;}
  const std::string omnesMatrixType() { return _omnesMatrixType;}
  
protected:
  //  virtual bool parseCommandLine(int argc,char **argv); 
  virtual bool parseCommandLine();
  std::string _configFile;
  std::string _keyName;
  unsigned int _noOfChannels;
  unsigned int _noOfPoles;
  std::vector<std::string> _gFactors;
  std::vector<std::string> _phpDescriptions;
  std::vector<std::string> _poles;
  std::string _projection;
  bool _useBarrierFactors;
  bool _useTruncatedBarrierFactors;
  unsigned int _orbitalMom;
  int _orderKMatrixBackground;
  int _orderPVectorBackground;
  bool _useAdler0;
  double _s0Adler;
  double _snormAdler;

  std::vector<std::string> _gFactorFixPoles;
  std::vector<std::string> _gFactorFixReleaseGFactor;
  std::vector<std::string> _gFactorFixKeepGRatio;
  std::vector<std::string> _gFactorFixMassPol;
  std::vector<std::string> _gFactorFixWidthPol;
  std::vector<std::string> _gFactorFixSeparateScale;

  std::map<std::string, double> _gFactorMap;
  std::map<std::string, double> _poleMap;

  std::map<std::pair<std::string, std::string>, std::string> _phpDescriptionMap;
  std::vector<std::pair<std::string, std::string>> _particlePairs;
  std::string _omnesMatrixType;
  
  po::options_description* _config;
};
