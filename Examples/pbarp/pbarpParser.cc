// Bertram Kopf (RUB)

#include "Examples/pbarp/pbarpParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


pbarpParser::pbarpParser(int argc,char **argv):
  ParserBase(argc,argv)
  ,_lMax(3)
  ,_pdgTableFile("/Particle/pdt.table")
{
  po::options_description common("Common Options");
  common.add_options()
    ("lmax", po::value<unsigned>(&_lMax)->default_value(_lMax),"choose lmax.")
    ("pdgTableFile",po::value<string>(&_pdgTableFile), "path of the pdg-table file relative to the top dir")
   ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ("finalStateParticle",po::value< vector<string> >(&_finalStateParticles),  "name of final state particles")
    ("production",po::value< vector<string> >(&_productionSystem),  "pair of produced particles")
    ("decay",po::value< vector<string> >(&_decaySystem),  "decay: mother and pair of decay particles")
    ("replaceParamSuffix",po::value< vector<string> >(&_replaceParSuffix),  "replace suffix for fit parameter name")
    ("histMass",po::value< vector<string> >(&_histMass),  "histograms inv mass for the selected final state paricles")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool pbarpParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);
  std::cout << "Maximum orbital momentum for pbarp system\t Lmax=\t" << _lMax <<std::endl;
  std::cout << "pdg table:\t" << _pdgTableFile <<std::endl;
  std::cout << "the final state particles are:" << std::endl;
  std::vector<std::string>::const_iterator it;
  for (it=_finalStateParticles.begin(); it!=_finalStateParticles.end();++it){
          std::cout << (*it) << "\n";
      }

  std::cout << "\nproduction system:" << std::endl;
  for (it=_productionSystem.begin(); it!=_productionSystem.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\ndecay system:" << std::endl;
  for (it=_decaySystem.begin(); it!=_decaySystem.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\nreplaced suffix for fit parameter name" << std::endl;
  for (it=_replaceParSuffix.begin(); it!=_replaceParSuffix.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\nhistograms inv mass for systems" << std::endl;
  for (it=_histMass.begin(); it!=_histMass.end();++it){
    std::cout << (*it) << "\n";
  }
  return true;
}
