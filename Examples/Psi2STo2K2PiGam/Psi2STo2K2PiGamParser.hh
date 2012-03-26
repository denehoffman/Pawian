/**
 * @file GArgumentParserLS.hpp
 */

/* Copyright (C) Dr. Ruediger Berlich and Karlsruhe Institute of Technology
 * (University of the State of Baden-Wuerttemberg and National Laboratory
 * of the Helmholtz Association)
 *
 * Contact: info [at] gemfony (dot) com
 *
 * This file is part of the Geneva library, Gemfony scientific's optimization
 * library.
 *
 * Geneva is free software: you can redistribute it and/or modify
 * it under the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */

// Standard headers go here

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Boost headers go here

#include <boost/version.hpp>

#if BOOST_VERSION < 103600
#error "Error: Boost should at least have version 1.36 !"
#endif /* BOOST_VERSION */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#ifndef Psi2STo2K2PiGamParser_HH
#define Psi2STo2K2PiGamParser_HH

// For Microsoft-compatible compilers
#if defined(_MSC_VER)  &&  (_MSC_VER >= 1020)
#pragma once
#endif

// GenEvA headers go here
//#include <common/GCommonEnums.hpp>
//#include <common/GSerializationHelperFunctionsT.hpp>
//#include <geneva/GOptimizationEnums.hpp>

namespace po = boost::program_options;
//using namespace Gem::Geneva;

class Psi2STo2K2PiGamParser
{
  public:
  typedef enum tagerrLogMode {debug,trace,routine,warning,error,alert} enErrLogMode;

  public:
    Psi2STo2K2PiGamParser(int argc,char **argv)
      : _configFile("./Psi2STo2K2PiGamConfig.cfg")
      , _errLogMode(debug)
      , _dataFile("/home/bertram/PawianGit110728/Pawian/Examples/Psi2STo2K2PiGam/data/110310_chic0_data_pwasample.dat")
      , _mcFile("/home/bertram/PawianGit110728/Pawian/Examples/Psi2STo2K2PiGam/data/110310_chic0_signalmc_pwasample.dat")
      , _hepMCinFile("/home/bertram/PawianGit110926/Pawian/Examples/Psi2STo2K2PiGam/HepMCEvt.in")
      , _paramFile("/Examples/Psi2STo2K2PiGam/data/startParamHyp8_2K0_2K2_K1_1270_woHyp5_KappaK1950V1.dat")
      , _startHypo("hyp8")
      , _qaMode(false)
      , _genMode(false)
      , _scanMode(false)
      , _branchingRatioMode(false)
    {
//     _disabledHyps.push_back("blainit");
    if (!parseCommandLine(argc, argv)) throw false;
  }
  
  const std::string& getConfigFile() const { return _configFile;}
  const enErrLogMode& getErrLogMode() const { return _errLogMode; }
  const std::string dataFile() const {return _dataFile;}
  const std::string mcFile() const {return _mcFile;}
  const std::string hepMCinFile() const {return _hepMCinFile;}
  const std::string fitParamFile() const {return _paramFile;}
  const std::vector<std::string>& enabledHyps() const { return _enabledHyps; }
  const std::string startHypo() const {return _startHypo;}
  const bool qaMode() const {return _qaMode;}
  const bool generatorMode() const {return _genMode;}
  const bool scanMode() const {return _scanMode;}
  const bool branchingRatioMode() const {return _branchingRatioMode;}
  const std::vector<std::string>& fixedParams() const { return _mnParFixs; }    
  const std::vector<std::string>& scanParams() const { return _mnParScan; }  
  const std::vector<std::string>& brParams() const { return _mnParCalcBR; }  
protected:
  bool parseCommandLine(int argc,char **argv);

protected:
  std::string _configFile;
  enErrLogMode _errLogMode;
  std::string _dataFile;
  std::string _mcFile;
  std::string _hepMCinFile;
  std::string _paramFile;
  std::string _startHypo;
  bool _qaMode;	
  bool _genMode;
  bool _scanMode;  
  bool _branchingRatioMode;
  std::vector<std::string> _enabledHyps;
  std::vector<std::string> _mnParFixs;
  std::vector<std::string> _mnParScan;
  std::vector<std::string> _mnParCalcBR;
};

#endif /* Psi2STo2K2PiGamParser_HH */
