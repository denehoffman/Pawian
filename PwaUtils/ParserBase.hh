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

class ParserBase {
public:
  typedef enum tagerrLogMode { debug, trace, routine, warning, error, alert } enErrLogMode;

  ParserBase(int argc,char **argv);
  virtual ~ParserBase(){;}
    
  const std::string& getConfigFile() const { return _configFile;}
  const enErrLogMode& getErrLogMode() const { return _errLogMode; }
  const std::string dataFile() const {return _dataFile;}
  const std::string mcFile() const {return _mcFile;}
  const std::string fitParamFile() const {return _paramFile;}
  const std::string outputFileNameSuffix() const {return _outputFileNameSuffix;}
  const std::vector<std::string>& enabledHyps() const { return _enabledHyps; }
  const std::string startHypo() const {return _startHypo;}
  const std::string mode() const {return _mode;}
  const std::vector<std::string>& fixedParams() const { return _mnParFixs; } 
  const int noOfThreads() const {return _noOfThreads;}
  const int ratioMcToData() const {return _ratioMcToData;}
  const bool cacheAmps() const {return _cacheAmps;}
  const bool useEvtWeight() const {return _useEvtWeight; }
  const std::string pdgTableFile() const {return _pdgTableFile;}
  
protected:
  virtual bool parseCommandLine(int argc,char **argv); 

  std::string _configFile;
  enErrLogMode _errLogMode;
  std::string _dataFile;
  std::string _mcFile;
  std::string _paramFile;
  std::string _startHypo;
  std::string _mode;		  
  std::string _outputFileNameSuffix;
  std::vector<std::string> _enabledHyps;
  std::vector<std::string> _mnParFixs;
  bool _verbose;
  int _noOfThreads;
  int _ratioMcToData;
  bool _cacheAmps;
  std::string _strErrLogMode;
  po::options_description* _desc;
  po::options_description* _common;
  po::options_description* _config;
  bool _useEvtWeight;
  std::string _pdgTableFile; 
};
