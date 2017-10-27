//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// ProdParamDynamics class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>
#include <math.h> 

#include "PwaUtils/ProdParamDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"

ProdParamDynamics::ProdParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::string type) :
  AbsDynamics(name, fsParticles, mother)
  ,_polOrder(0)
{
  if(type=="FormPol0") _polOrder=0;
  else if(type=="FormPol1") _polOrder=1;
  else if(type=="FormPol2") _polOrder=2;
  else{
    Alert <<"dynamics of the type" << name << "is not supported!!!" << endmsg;
    exit(1); 
  }
  _currentPolParams.resize(_polOrder+1);
  _fitPolParNames.resize(_polOrder+1);
  _isLdependent=false;
}

ProdParamDynamics::~ProdParamDynamics()
{
}

complex<double> ProdParamDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  //  complex<double> result=BreitWignerFunction::NonRel(theData->DoubleMassId.at(_dynId), _currentMass, _currentWidth); 
  complex<double> result(0.,0.);
  for(unsigned int i=0; i<=_polOrder; ++i){
    complex<double> currentResult(_currentPolParams.at(i)*pow(theData->DoubleMassId.at(_dynId), i), 0.);
    result+=currentResult;
  }
  result*=exp(_currentExpParam*theData->DoubleMassId.at(_dynId));
 
  if ( _cacheAmps){
     theMutex.lock();
     _cachedMap[evtNo]=result;
     theMutex.unlock();
  } 
  

  return result;
}


void ProdParamDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  
  for(unsigned int i=0; i<=_polOrder; ++i){
    std::stringstream stringStr;
    stringStr << i;
    std::string fitParName="ProdPol"+stringStr.str();
    _fitPolParNames[i]=fitParName;
    if(i==0){
      fitPar->Add(_fitPolParNames.at(i), 1., 0.1);
      _currentPolParams[i]=1.;
    }
    else{
      fitPar->Add(_fitPolParNames.at(i), 0., 0.1);
      _currentPolParams[i]=0.;
    } 
  }
 fitPar->Add("ExpProd", 0., 0.1);
 _currentExpParam=0.;
}


void ProdParamDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  for(unsigned int i=0; i<=_polOrder; ++i){
    _currentPolParams.at(i)=fitPar->Value(_fitPolParNames.at(i));
  }
  _currentExpParam=fitPar->Value("ExpProd");
}

void ProdParamDynamics::fillParamNameList(){
  _paramNameList.clear();
  for(unsigned int i=0; i<=_polOrder; ++i){
    std::stringstream stringStr;
    stringStr << i;
    std::string fitParName="ProdPol"+stringStr.str();
    _fitPolParNames[i]=fitParName;
    _paramNameList.push_back(fitParName);
  }
  _paramNameList.push_back("ExpProd");
}

// bool ProdParamDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
//   std::vector<std::string>::iterator it;
//   for(it=_paramNameList.begin(); it!=_paramNameList.end(); ++it){
//     if(
//   } 
// }
