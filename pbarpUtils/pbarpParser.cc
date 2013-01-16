// pbarpParser class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf


#include "pbarpUtils/pbarpParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


pbarpParser::pbarpParser(int argc,char **argv):
  ParserBase(argc,argv)
  ,_lMax(3)
{
  po::options_description common("Common Options");
  common.add_options()
    ("lmax", po::value<unsigned>(&_lMax)->default_value(_lMax),"choose lmax.")
    ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    //    ("finalStateParticle",po::value< vector<string> >(&_finalStateParticles),  "name of final state particles")
    ("production",po::value< vector<string> >(&_productionSystem),  "pair of produced particles")
    ("histMass",po::value< vector<string> >(&_histMass),  "histograms inv mass for the selected final state paricles")
    ("histAngles",po::value< vector<string> >(&_histAngles),  "histograms decay angles")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool pbarpParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);
  std::cout << "Maximum orbital momentum for pbarp system\t Lmax=\t" << _lMax <<std::endl;

  std::vector<std::string>::const_iterator it;
  std::cout << "\nproduction system:" << std::endl;
  for (it=_productionSystem.begin(); it!=_productionSystem.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\nhistograms inv mass for systems" << std::endl;
  for (it=_histMass.begin(); it!=_histMass.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\nhistograms decay angles for systems" << std::endl;
  for (it=_histAngles.begin(); it!=_histAngles.end();++it){
    std::cout << (*it) << "\n";
  }
  return true;
}
