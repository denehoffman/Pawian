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

// Bertram Kopf (RUB)

#include "ConfigParser/KMatrixParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include <iterator>
#include <iostream>
#include <fstream>
using namespace std;


KMatrixParser::KMatrixParser(std::string& path)
      : _configFile(path)
      , _keyName("Default")
      , _noOfChannels(0)
      , _noOfPoles(0)
      , _projection("")
      , _useBarrierFactors(false)
      ,_useTruncatedBarrierFactors(false)
      , _orbitalMom(0)
      ,_orderKMatrixBackground(-1)
      ,_orderPVectorBackground(-1)
      ,_useAdler0(false)
      ,_s0Adler(1.)
      ,_snormAdler(1.)
      ,_omnesMatrixType("PiPiS5Channel")
      ,_fixedParamFile("")
      ,_config(new po::options_description("Configuration file options"))
     {
    _config->add_options()
      ("keyName",po::value<std::string>(&_keyName), "name for this specific k-matrix")
      ("noOfChannels",po::value<unsigned int>(&_noOfChannels),  "number of channels")
      ("noOfPoles",po::value<unsigned int>(&_noOfPoles),  "number of poles")
      ("gFactor",po::value< vector<string> >(&_gFactors),  "g-factors for the individual channels")
      ("phpDescription",po::value< vector<string> >(&_phpDescriptions),  "description of the phasespace factors for the individual channels") 
      ("pole",po::value< vector<string> >(&_poles),  "pole: name and mass")
      ("projection",po::value< string>(&_projection),  "projection of the P-vector via pair of decay particles")
      ("useBarrierFactors",po::value<bool>(&_useBarrierFactors),  "calculation with or without barrier factors")
      ("useTruncatedBarrierFactors",po::value<bool>(&_useTruncatedBarrierFactors), "use truncated barrier factors, needed for tensor formalism")
      ("orbitalMomentum",po::value<unsigned int>(&_orbitalMom),  "orbital momentum of the decay")
      ("orderKMatrixBackground",po::value<int>(&_orderKMatrixBackground),  "order of the K-matrix background")
      ("orderPVectorBackground",po::value<int>(&_orderPVectorBackground),  "order of the P-vector background")
      ("useAdler0",po::value<bool>(&_useAdler0),  "use adler0 term")
      ("s0Adler",po::value<double>(&_s0Adler),  "s0Adler parameter")
      ("snormAdler",po::value<double>(&_snormAdler),  "snormAdler parameter")
      ("gFactorFixPoles",po::value< vector<string> >(&_gFactorFixPoles), "poles with fixed g-factor ratios, total widths and total mass")
      ("gFactorFixReleaseGFactor",po::value< vector<string> >(&_gFactorFixReleaseGFactor), "free g-factor for fixing g-factor ratios, total widths and total mass")
      ("gFactorFixKeepGRatio",po::value< vector<string> >(&_gFactorFixKeepGRatio), "g-factors to be scaled for fixing g-factor ratios, total widths and total mass")
      ("gFactorFixMassPol",po::value< vector<string> >(&_gFactorFixMassPol), "3rd order polynomial to correct the pole mass for fixing g-factor ratios, total widths and total mass")  
      ("gFactorFixWidthPol",po::value< vector<string> >(&_gFactorFixWidthPol), "3rd order polynomial to correct the pole width for fixing g-factor ratios, total widths and total mass")
      ("gFactorFixSeparateScale",po::value< vector<string> >(&_gFactorFixSeparateScale), "g-factors to be scaled separately related to free g-factor")
      ("omnesMatrixType",po::value<string>(&_omnesMatrixType),"type of the Omnes-matrix")
      ("fixedParamFile",po::value<string>(&_fixedParamFile), "file with K-Matrix parameters to be used for a fixed parametrization")
      ;
    parseCommandLine();
  }

/************************************************************************************************/
/************************************************************************************************/
/**
  * A function that parses the command line for all required parameters
  */

bool KMatrixParser::parseCommandLine()
{

    std::ifstream ifs(_configFile.c_str());
    if(!ifs.good())
      {
	std::cout << "Error accessing configuration file " << _configFile << std::endl;
 	exit(0);
      }

     po::options_description config_file_options;
     config_file_options.add(*_config);

    po::variables_map vm;
    store(po::parse_config_file(ifs, config_file_options), vm);
    po::notify(vm);

    std::cout << "\nRunning with the following options using " << _configFile << ":\n\n"
	      << "key name: " << _keyName  << "\n\n"
	      << "number of channels: " << _noOfChannels  << "\n\n"
	      << "number of poles: " << _noOfPoles  << "\n\n"
	      << "\nprojection: " << _projection << "\n\n"
	      << "with barrier factors: " << _useBarrierFactors  << "\n\n"
	      << "truncated barrier factors: " << _useTruncatedBarrierFactors  << "\n\n"
	      << "orbital momentum: " << _orbitalMom  << "\n\n"
	      << "order of the K-matrix background: " << _orderKMatrixBackground << "\n\n"
	      << "order of the P-vector background: " << _orderPVectorBackground << "\n\n"
	      << "use Adler0 term: " << _useAdler0 << "\n\n"
	      << "s0Adler: " << _s0Adler << "\n\n"
	      << "snormAdler: " << _snormAdler << "\n\n"
	      << "fixedParamFile: " << _fixedParamFile  << "\n\n"
	      << endl;

    std::cout << "g-factors are defined as follows:" << "size: " << _gFactors.size() << std::endl;
    std::vector<std::string>::const_iterator it;
    for(it=_gFactors.begin(); it!=_gFactors.end(); ++it){
      std::cout << (*it) << "\n";
      std::istringstream particles(*it);
      std::string firstParticleName;
      std::string secondParticleName;
      particles >> firstParticleName >> secondParticleName;
      std::pair<std::string, std::string > currentParticlePair=make_pair(firstParticleName,secondParticleName);
      _particlePairs.push_back(currentParticlePair);
      std::vector<std::string>::const_iterator itPoles;
      for(itPoles=_poles.begin(); itPoles!=_poles.end(); ++itPoles){

	std::istringstream poleIString(*itPoles);
	std::string currentPoleName;

	poleIString >> currentPoleName;

	std::string currentgFacName=currentPoleName+firstParticleName+secondParticleName+"gFactor";

	std::string currentgValueStr;
	if(!(particles >> currentgValueStr)){
	  Alert << "g-factors " << currentgFacName << " does not exist!" << endmsg;
	  exit(0);
	}
	std::istringstream currentgValueiStr(currentgValueStr);
	double currentGValue;
	if (!(currentgValueiStr >> currentGValue)){
	  Alert << "cannot convert " << currentgValueStr << " to a double value" << endmsg;
	  exit(0);
	}

	_gFactorMap[currentgFacName]=currentGValue;
      }
    }


    for(it=_phpDescriptions.begin(); it!=_phpDescriptions.end(); ++it){
      std::istringstream particles(*it);
      std::string firstParticleName;
      std::string secondParticleName;
      std::string phpDescription;
      particles >> firstParticleName >> secondParticleName >> phpDescription;
      std::pair<std::string, std::string > particlePair=make_pair(firstParticleName,secondParticleName);
      _phpDescriptionMap[particlePair]=phpDescription;        
      }

   // fill remaining parts with default description
   if ( _phpDescriptions.size() < _particlePairs.size() ){ //fill remaining with default description
     std::vector< std::pair<std::string, std::string> >::iterator itPPair;
     for(itPPair= _particlePairs.begin(); itPPair!= _particlePairs.end(); ++itPPair){
     // check if phpDescription already exists
     std::map<std::pair<std::string, std::string>, std::string>::iterator found=_phpDescriptionMap.find(*itPPair);
     if ( !(found!=_phpDescriptionMap.end()) ) _phpDescriptionMap[*itPPair]="Default";
   }
}
  

   std::cout << "\nthe descriptions for the phase space factors are defined as follows:" << std::endl;
   std::map<std::pair<std::string, std::string>, std::string>::iterator itPhpDesc;
   for(itPhpDesc= _phpDescriptionMap.begin(); itPhpDesc!= _phpDescriptionMap.end(); ++itPhpDesc){
   std::cout << itPhpDesc->first.first << "\t" << itPhpDesc->first.second << "\t" << itPhpDesc->second << std::endl; 
}

    std::cout << "\npoles: name and mass" << std::endl;
    for(it=_poles.begin(); it!=_poles.end(); ++it){
      std::cout << (*it) << "\n";
    }

    std::map<std::string, double>::iterator itgMap;
    for (itgMap=_gFactorMap.begin(); itgMap!=_gFactorMap.end(); ++itgMap){
      std::cout << "\n"<< itgMap->first  << "\t" << itgMap->second;      
    }
    std::cout << "\n"; 
  return true;
}

