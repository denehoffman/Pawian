//************************************************************************//
//									  //
//  Copyright 2018 Bertram Kopf (bertram@ep1.rub.de)			  //
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

#include <getopt.h>
#include <fstream>
#include <sstream>
#include <string>
//#include <boost/multi_array.hpp>
#include "KMatrixExtract/TMatrixExtrBase.hh"
#include "qft++/topincludes/relativistic-quantum-mechanics.hh" 
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/TMatrixBase.hh"
#include "PwaDynamics/TMatrixRel.hh"
#include "PwaDynamics/TMatrixNonRel.hh"
#include "PwaDynamics/KMatrixBase.hh"
#include "PwaDynamics/KPole.hh"
#include "PwaDynamics/KPoleBarrier.hh"
#include "PwaDynamics/KMatrixRel.hh"
#include "PwaDynamics/KMatrixRelBg.hh"
#include "PwaDynamics/AbsPhaseSpace.hh"
#include "PwaDynamics/PhaseSpaceFactory.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/PawianConstants.hh"

#include "ConfigParser/ParserBase.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "FitParams/ParamFactory.hh"

#include "ErrLogger/ErrLogger.hh"

TMatrixExtrBase::TMatrixExtrBase(std::string pathToConfigParser, std::string pathToFitParams, std::string sheet) :
  _kMatrixParser(new KMatrixParser(pathToConfigParser))
  ,_pathToFitParams(pathToFitParams)
  ,_sheet(sheet)
  ,_orbitalL(0)
{ 
  init();
}

TMatrixExtrBase::~TMatrixExtrBase()
{
}

void TMatrixExtrBase::init(){
  GlobalEnv::instance()->setup();
  _particleTable=GlobalEnv::instance()->particleTable();
  _tMatDynPtr= std::shared_ptr<TMatrixDynamics>(new TMatrixDynamics(_kMatrixParser));

    std::shared_ptr<AbsPawianParameters> params=ParamFactory::instance()->getParametersPointer("Pawian");
	_tMatDynPtr->fillDefaultParams(params);
	_tMatDynPtr->fillParamNameList();
	_kMatrixParamNames =_tMatDynPtr->paramNames();
	InfoMsg << "kMatrixParamNames.size(): " << _kMatrixParamNames.size() << endmsg;

	std::ifstream ifs(_pathToFitParams);
	if(!ifs.good()) 
	{ //file doesn't exist; dum default params
	  WarningMsg << "could not parse " << _pathToFitParams << endmsg;
	  WarningMsg << "dump defaut parameter " << _pathToFitParams << endmsg;
	  std::string defaultparamsname="defaultParams.dat";
	  std::ofstream theStreamDefault ( defaultparamsname );
	  params->print(theStreamDefault);
	  theStreamDefault.close();
	  exit(1);        
	}   


	AbsPawianParamStreamer thePawianStreamer(_pathToFitParams);
	_params = thePawianStreamer.paramList();

	InfoMsg << "The k-Matrix input parameter are: " << endmsg;
	_params->print(std::cout);
	if(_pathToFitParams != "") _tMatDynPtr->updateFitParams(_params); 

	_kMatr = _tMatDynPtr->getKMatix(); 
	_tMatr = _tMatDynPtr->getTMatix();  

	_phpVecs=_kMatr->phaseSpaceVec();
	_gFactorNames= _tMatDynPtr->gFactorNames();
	_orbitalL= _tMatDynPtr->orbitalL();

	std::vector<std::string> poleNameAndMassVecs=_kMatrixParser->poles();
	std::vector<std::string>::iterator itString;
	for (itString=poleNameAndMassVecs.begin(); itString!=poleNameAndMassVecs.end(); ++itString){
	  std::istringstream poleIString(*itString);
	  std::string currentPoleName;
	  std::string currentPoleMassStr;
	  poleIString >> currentPoleName >> currentPoleMassStr;

	  std::istringstream currentPoleMassiStr(currentPoleMassStr);
	  double currentValue;
	  if(!(currentPoleMassiStr >> currentValue)){
		Alert << "cannot convert " << currentPoleMassStr << " to a double value" << endmsg;
		exit(0);
	  }
	}

	const std::vector<std::string> gFacStringVec=_kMatrixParser->gFactors();
	DebugMsg << "gFacStringVec.size(): " << gFacStringVec.size() << endmsg;
	_signs.resize(gFacStringVec.size());
	for (int idx=0; idx<gFacStringVec.size(); ++idx) _signs[idx]=-1.;
	unsigned int iSize = _sheet.size();
	for(unsigned int i = 0; i < iSize; i++)
	{
	  if(_sheet[i]=='p')  {  _signs[i] = 1;} 
	}
	std::map<std::pair<std::string, std::string>, std::string> phpDescriptionVec=_kMatrixParser->phpDescriptionMap();
	std::cout << "phpDescriptionVec.size(): " << phpDescriptionVec.size() << std::endl;
	std::vector<std::string>::const_iterator itStrConst;
	for(itStrConst=gFacStringVec.begin(); itStrConst!=gFacStringVec.end(); ++itStrConst){
	  std::istringstream particles(*itStrConst);
	  std::string firstParticleName;
	  std::string secondParticleName;
	  particles >> firstParticleName >> secondParticleName;
	  std::pair<std::string, std::string> currentParticlePair=make_pair(firstParticleName, secondParticleName);
	  Particle* firstParticle = _particleTable->particle(firstParticleName);
	  Particle* secondParticle = _particleTable->particle(secondParticleName);
	  if(0==firstParticle || 0==secondParticle){
		Alert << "particle with name: " << firstParticleName <<" or " << secondParticleName << " doesn't exist in pdg-table" << endmsg;
		exit(0);
	  }
	}
	_tMatr->SetBumImPartSigns(_signs);
	return;
}

void TMatrixExtrBase::updateTMatDy(std::shared_ptr<AbsPawianParameters> params) {
  _tMatDynPtr->updateFitParams(params); 

  _kMatr = _tMatDynPtr->getKMatix(); 
  _tMatr = _tMatDynPtr->getTMatix();  
  return;
}
