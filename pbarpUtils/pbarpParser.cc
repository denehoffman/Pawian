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
  ,_pbarMomentum(2.)
{
  po::options_description common("Common Options");
  common.add_options()
    ("lmax", po::value<unsigned>(&_lMax)->default_value(_lMax),"choose lmax.")
    ("pbarmom", po::value<float>(&_pbarMomentum)->default_value(_pbarMomentum),"antiproton momentum")
    ;
  
  _common->add(common);
  
  po::options_description config("Configuration file options");
  config.add_options()
    ("histMass",po::value< vector<string> >(&_histMass),  "histograms inv mass for the selected final state paricles")
    ("histAngles",po::value< vector<string> >(&_histAngles),  "histograms decay angles")
    ("spinDensity", po::value< vector<string> >(&_spinDensity), "particles for spin density matrix calculation")
    ("lmaxParticle", po::value< vector<string> >(&_lmaxParticle), "individual lmax (pbarp system) for certain particle")
    ;
  _config->add(config);
  
  parseCommandLine(argc, argv);
}

bool pbarpParser::parseCommandLine(int argc, char **argv)
{
  ParserBase::parseCommandLine(argc, argv);

  std::cout << "Antiproton momentum in lab frame\t pbarmom = " << _pbarMomentum << std::endl;
  std::cout << "Maximum orbital momentum for pbarp system\t Lmax = " << _lMax <<std::endl;

  std::vector<std::string>::const_iterator it;
  // std::cout << "\nproduction system:" << std::endl;
  // for (it=_productionSystem.begin(); it!=_productionSystem.end();++it){
  //   std::cout << (*it) << "\n";
  // }

  std::cout << "\nhistograms inv mass for systems" << std::endl;
  for (it=_histMass.begin(); it!=_histMass.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\nhistograms decay angles for systems" << std::endl;
  for (it=_histAngles.begin(); it!=_histAngles.end();++it){
    std::cout << (*it) << "\n";
  }

  std::cout << "\nspin density matrix calculation for particles" << std::endl;
  for (it=_spinDensity.begin(); it!=_spinDensity.end(); ++it){
     std::cout << (*it) << "\n";
  }

  std::cout << "\nindividual lmax settings " << std::endl;
  for (it=_lmaxParticle.begin(); it!=_lmaxParticle.end(); ++it){
     std::cout << (*it) << "\n";
  }

  std::cout << std::endl;

  return true;
}
