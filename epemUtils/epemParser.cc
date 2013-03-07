// epemParser class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf


#include "epemUtils/epemParser.hh"
#include "ErrLogger/ErrLogger.hh"

#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


epemParser::epemParser(int argc,char **argv):
  ParserBase(argc,argv)
{
  po::options_description common("Common Options");
  common.add_options()
     ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ("spinDensity", po::value< vector<string> >(&_spinDensity), "particles for spin density matrix calculation") 
    ("cmsMass", po::value<double>(&_cmsMass), "CMS mass")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool epemParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);

  std::cout << "\ncms mass:\t" << _cmsMass << std::endl;

  std::vector<std::string>::const_iterator it;
  
  std::cout << "\nspin density matrix calculation for particles" << std::endl;
   for (it=_spinDensity.begin(); it!=_spinDensity.end(); ++it){
      std::cout << (*it) << "\n";
   }
  
  std::cout << std::endl;
  
  return true;
}
