#ifndef JpsiGamEtaPiPiParser_HH
#define JpsiGamEtaPiPiParser_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

// Boost headers go here

#include <boost/version.hpp>

#if BOOST_VERSION < 103600
#error "Error: Boost should at least have version 1.36 !"
#endif /* BOOST_VERSION */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

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

class JpsiGamEtaPiPiParser
{
  public:
  typedef enum tagerrLogMode {debug,trace,routine,warning,error,alert} enErrLogMode;

  public:
    JpsiGamEtaPiPiParser(int argc,char **argv)
      : _configFile("./JpsiGamEtaPiPi.cfg")
      , _errLogMode(debug)
      , _dataFile("")
      , _mcFile("")
      , _paramFile("")
      , _startHypo("base")
      , _mode("plotmode")
      , _massIndependentFit(false)
      , _useCommonProductionPhases(false)
      , _massMin(0.7)
      , _massMax(3.1)
  {
    if (!parseCommandLine(argc, argv)) throw false;
  }
  
  const std::string& getConfigFile() const { return _configFile;}
  const enErrLogMode& getErrLogMode() const { return _errLogMode; }
  const std::string dataFile() const {return _dataFile;}
  const std::string mcFile() const {return _mcFile;}
  const std::string fitParamFile() const {return _paramFile;}
  const std::vector<std::string>& enabledHyps() const { return _enabledHyps; }
  const std::string startHypo() const {return _startHypo;}
  const std::string mode() const {return _mode;}
  const std::vector<std::string>& fixedParams() const { return _mnParFixs; } 
  const bool massIndependentFit() const {return _massIndependentFit; }
  const bool useCommonProductionPhases() const {return _useCommonProductionPhases; }
  const std::pair<double, double> massRange() const { return std::make_pair( _massMin, _massMax  ) ; }
  
  
protected:
  bool parseCommandLine(int argc,char **argv);

protected:
  std::string _configFile;
  enErrLogMode _errLogMode;
  std::string _dataFile;
  std::string _mcFile;
  std::string _paramFile;
  std::string _startHypo;
  std::string _mode;		  
  std::vector<std::string> _enabledHyps;
  std::vector<std::string> _mnParFixs;
  bool _massIndependentFit;
  bool _useCommonProductionPhases;
  
  double _massMin;
  double _massMax;

};


#endif /* JpsiGamEtaPiPiParser_HH */
