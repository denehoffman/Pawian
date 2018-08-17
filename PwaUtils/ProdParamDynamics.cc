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
#include "PwaUtils/GlobalEnv.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"
#include "Utils/IdStringMapRegistry.hh"

ProdParamDynamics::ProdParamDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, ChannelID channelID, std::string type) :
  AbsDynamics(name, fsParticles, mother)
  ,_polOrder(0)
  ,_channelID(channelID)
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
  //  int evtNo=theData->evtNo;
  if (!_recalculate){
    return _cachedMap[theData->evtNo];
  }

  //  complex<double> result=BreitWignerFunction::NonRel(theData->DoubleMassId.at(_dynId), _currentMass, _currentWidth); 
  complex<double> result(0.,0.);
  for(unsigned int i=0; i<=_polOrder; ++i){
    complex<double> currentResult(_currentPolParams.at(i)*pow(theData->DoubleMassId.at(_dynId), i), 0.);
      // InfoMsg << "pol(" << i << "): " << currentResult
      // 	      <<"\ttheData->DoubleMassId.at(" << _dynId << "): " << theData->DoubleMassId.at(_dynId) 
      // 	      << "\tcurrentPolParams.at(i): " << _currentPolParams.at(i) 
      // 	      << "\tpow(theData->DoubleMassId.at(_dynId), i), 0.): " << pow(theData->DoubleMassId.at(_dynId), i) 
      // 	      << endmsg; 
      

    result+=currentResult;
  }
  result*=exp(_currentExpParam*theData->DoubleMassId.at(_dynId));
 
  if ( _cacheAmps){
     theMutex.lock();
     _cachedMap[theData->evtNo]=result;
     theMutex.unlock();
  } 
  

  return result;
}


void ProdParamDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  fillParamNameList();  
  for(unsigned int i=0; i<=_polOrder; ++i){
    //   std::stringstream stringStr;
    //    stringStr << i;
    //   std::string fitParName="ProdPol"+stringStr.str();
    //    _fitPolParNames[i]=fitParName;
    if(i==0){
      fitPar->Add(_fitPolParNames.at(i), 1., 0.1);
      _currentPolParams[i]=1.;
    }
    else{
      fitPar->Add(_fitPolParNames.at(i), 0., 0.1);
      _currentPolParams[i]=0.;
    } 
  }
  // fitPar->Add("ExpProd", 0., 0.1);
 fitPar->Add(_fitExpParName, 0., 0.1);
 _currentExpParam=0.;
}


void ProdParamDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  for(unsigned int i=0; i<=_polOrder; ++i){
    _currentPolParams.at(i)=fitPar->Value(_fitPolParNames.at(i));
//    InfoMsg << "updated par " << _fitPolParNames.at(i) << ": " << fitPar->Value(_fitPolParNames.at(i)) << endmsg;
  }
  _currentExpParam=fitPar->Value(_fitExpParName);
//  InfoMsg << "updated par " << _fitExpParName << ": " << fitPar->Value(_fitExpParName) << endmsg;
}

void ProdParamDynamics::fillParamNameList(){
  _paramNameList.clear();
  std::stringstream stringStrChannelId;
    stringStrChannelId << _channelID;
  for(unsigned int i=0; i<=_polOrder; ++i){
    std::stringstream stringStr;
    stringStr << i;
    std::string fitParName="ProdPol"+stringStr.str()+"Channel"+stringStrChannelId.str();
    _fitPolParNames[i]=fitParName;
    _paramNameList.push_back(fitParName);
  }
  _fitExpParName="ExpProdChannel"+stringStrChannelId.str();
  _paramNameList.push_back(_fitExpParName);
}

// bool ProdParamDynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
//   std::vector<std::string>::iterator it;
//   for(it=_paramNameList.begin(); it!=_paramNameList.end(); ++it){
//     if(
//   } 
// }

void ProdParamDynamics::fillMasses(EvtData* theData){
  std::vector<Particle*> theFinalStateParticles=GlobalEnv::instance()->Channel(_channelID)->finalStateParticles();
  Vector4<double> mass4Vec(0.,0.,0.,0.);
  std::vector<Particle*>::iterator it;
  for (it=theFinalStateParticles.begin(); it != theFinalStateParticles.end(); ++it){
    std::string currentName=(*it)->name();
    mass4Vec+=theData->FourVecsId.at(IdStringMapRegistry::instance()->stringId(currentName));
  }

  double sqrMass=mass4Vec.Mass2();
  if (sqrMass>0.) theData->DoubleMassId[_dynId]=mass4Vec.Mass();
  else if( sqrMass > -1.e-6) theData->DoubleMassId[_dynId]=0.;
  else{
    Alert << "mass4Vec.Mass2() is " << mass4Vec.Mass2() << " and thus < -1e-6 !!!" 
	  << "\nexit !!!" << endmsg;
    exit(0); 
  }
}
