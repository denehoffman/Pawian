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
#include "Particle/PdtParser.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "FitParams/AbsPawianParamStreamer.hh"
#include "FitParams/AbsPawianParameters.hh"
#include "Utils/PawianConstants.hh"
#include "ConfigParser/ParserBase.hh"
#include "ConfigParser/KMatrixParser.hh"
#include "ConfigParser/pipiScatteringParser.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/TMatrixDynamics.hh"
#include "FitParams/ParamFactory.hh"
#include "ErrLogger/ErrLogger.hh"

TMatrixExtrBase::TMatrixExtrBase(pipiScatteringParser* theParser) :
  TMatrixGeneralBase(theParser)
  ,_sheet(theParser->sheet()){

  InfoMsg << "_sheet.size(): " << _sheet.size() << endmsg;
  _signs.resize(_sheet.size());
  _sheetSize=_signs.size();
  setSheet(_sheet);
  fillParams();
}

TMatrixExtrBase::~TMatrixExtrBase()
{
}

void TMatrixExtrBase::setSheet(std::string newSheet){
  if(newSheet.size() != _signs.size()){
    ErrMsg << "size of the new sheet: " << newSheet.size() << " is not in agreement with the size of the T-Matrix: " <<  _signs.size() << endmsg;
    exit(1);
  }

  for(unsigned int i = 0; i < _signs.size(); i++) {
    if(newSheet[i]=='p' || newSheet[i]=='+')  {  _signs[i] = 1;}
    else if(newSheet[i]=='n' || newSheet[i]=='-') _signs[i] = -1;
    else{
      ErrMsg << "sheet with name " << newSheet[i] << " is not supported!!!" << endmsg;
      exit(1);
    }
  }
  _tMatr->SetBumImPartSigns(_signs);
}

double TMatrixExtrBase::calcTMatrix(double eReal, double eImag){
  _tMatr->evalMatrix( std::complex<double>(eReal, eImag), _orbitalL );
  return std::abs((*_tMatr)(_decProjectionIndex,_decProjectionIndex));
}

void TMatrixExtrBase::updateTMatDy(std::shared_ptr<AbsPawianParameters> params) {
  _tMatrDyn->updateFitParams(params);
  _kMatr = _tMatrDyn->getKMatix(); 
  _tMatr = _tMatrDyn->getTMatix();
  //  _tMatr->SetBumImPartSigns(_signs);  
  return;
}

std::shared_ptr<TMatrixRel> TMatrixExtrBase::getNewTMat(){
 std::shared_ptr<TMatrixDynamics> currentTMatDyn(new TMatrixDynamics(_kMatrixParser));
 currentTMatDyn->updateFitParams(_params);
 std::shared_ptr<TMatrixRel> currentTMat=currentTMatDyn->getTMatix();
 currentTMat->SetBumImPartSigns(_signs);
 return currentTMat; 
}

void TMatrixExtrBase::SetParamValue(const std::string & paramName, double paramVal){
  TMatrixGeneralBase::SetParamValue(paramName, paramVal);
  //  _params->print(std::cout);
  updateTMatDy(_params);
}

int TMatrixExtrBase::sheetSize(){
  return _signs.size();
}
