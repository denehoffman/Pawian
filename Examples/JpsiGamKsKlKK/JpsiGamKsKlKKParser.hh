#ifndef JpsiGamKsKlKKParser_HH
#define JpsiGamKsKlKKParser_HH

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

class JpsiGamKsKlKKParser
{
  public:
  typedef enum tagerrLogMode {debug,trace,routine,warning,error,alert} enErrLogMode;

  public:
    JpsiGamKsKlKKParser(int argc,char **argv)
      : _configFile("./JpsiGamKsKlKKConfig.cfg")
      , _errLogMode(debug)
      , _dataFile("/data/sleipnir1/bertram/JpsiGamKsKlKKData/fitvectorDATA_kkpi0.dat")
      , _mcFile("/data/sleipnir1/bertram/JpsiGamKsKlKKData/fitvectorMC_chic1_kkpi0.dat")
      , _paramFile("/data/sleipnir1/bertram/JpsiGamKsKlKKData/startParamSpin02.dat")
      , _startHypo("base")
      , _qaMode(false)
    {
    if (!parseCommandLine(argc, argv)) throw false;
  }
  
  const std::string& getConfigFile() const { return _configFile;}
  const enErrLogMode& getErrLogMode() const { return _errLogMode; }
  const std::string dataFile() const {return _dataFile;}
  const std::string mcFile() const {return _mcFile;}
  const std::string fitParamFile() const {return _paramFile;}
  const std::vector<std::string>& disabledHyps() const { return _disabledHyps; }
  const std::string startHypo() const {return _startHypo;}
  const bool qaMode() const {return _qaMode;}
  const std::vector<std::string>& fixedParams() const { return _mnParFixs; }    
protected:
  bool parseCommandLine(int argc,char **argv);

protected:
  std::string _configFile;
  enErrLogMode _errLogMode;
  std::string _dataFile;
  std::string _mcFile;
  std::string _paramFile;
  std::string _startHypo;
  bool _qaMode;		  
  std::vector<std::string> _disabledHyps;
  std::vector<std::string> _mnParFixs;
};


#endif /* JpsiGamKsKlKKParser_HH */
