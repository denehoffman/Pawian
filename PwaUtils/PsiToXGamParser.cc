// Bertram Kopf (RUB)

#include "PwaUtils/PsiToXGamParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


PsiToXGamParser::PsiToXGamParser(int argc,char **argv):
  ParserBase(argc,argv)
  ,_massIndependentFit(false)
  , _useCommonProductionPhases(false)
  , _massMin(0.7)
  , _massMax(3.1)
  , _jobOption("")
{
  po::options_description common("Common Options");
  common.add_options()
    ("massIndependentFit", po::value<bool>(&_massIndependentFit), "enable/disable mass independence in fit")
    ("commonProdPhases",po::value<bool>(&_useCommonProductionPhases), "enable/disable common production phases")
    ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ("massRangeMin",po::value<double>(&_massMin), "min of eta pi pi mass range for mass indep. fit")
    ("massRangeMax",po::value<double>(&_massMax), "max of eta pi pi mass range for mass indep. fit")
    ("jobOption",po::value<string>(&_jobOption), "add a jobOption string to all produced files")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

