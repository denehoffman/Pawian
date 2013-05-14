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
     ("spinDensity", po::value< vector<string> >(&_spinDensity), "particles for spin density matrix calculation")
     ("lmaxParticle", po::value< vector<string> >(&_lmaxParticle), "individual lmax (pbarp system) for certain particle")
     ("dropPbarpLForParticle", po::value< vector<string> >(&_dropPbarpLForParticle), "drop l (pbarp system) for certain particle")
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

  std::cout << "\nspin density matrix calculation for particles" << std::endl;
  for (it=_spinDensity.begin(); it!=_spinDensity.end(); ++it){
     std::cout << (*it) << "\n";
  }

  std::cout << "\nindividual pbarp angular momentum settings " << std::endl;
  for (it=_lmaxParticle.begin(); it!=_lmaxParticle.end(); ++it){
     std::cout << (*it) << "\n";
  }

  for (it=_dropPbarpLForParticle.begin(); it!=_dropPbarpLForParticle.end(); ++it){
     std::cout << (*it) << "\n";
  }

  std::cout << std::endl;

  return true;
}
