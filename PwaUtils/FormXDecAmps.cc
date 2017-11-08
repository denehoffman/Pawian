//************************************************************************//
//									  //
//  Copyright 2017 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// FormXDecAmps class definition file. -*- C++ -*-
// Copyright 2017 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/FormXDecAmps.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/FormationDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

FormXDecAmps::FormXDecAmps(std::shared_ptr<FormationDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
{
  initialize();
}

FormXDecAmps::FormXDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
{
  initialize();
}

FormXDecAmps::~FormXDecAmps()
{
}

complex<double> FormXDecAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){
  //dummy
  complex<double> result(0.,0.);
  return result;
}


complex<double> FormXDecAmps::XdecAmp(Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp){
  complex<double> result(0.,0.);
  if(_J<fabs(lamX)) return result;

  int evtNo=theData->evtNo;

  short currentSpinIndex=FunctionUtils::spin1IdIndex(_projIdThreadMap.at(std::this_thread::get_id()),lamX);
  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpIdMap.at(evtNo).at(_absDyn->grandMaId(grandmaAmp)).at(currentSpinIndex);
    return result;
  }

  Spin absLamX(lamX);
  if(lamX<0) absLamX=-lamX;
  theMutex.lock();
  complex<double> currentDyn=_cachedDynIdMap.at(std::this_thread::get_id()).at(_absDyn->grandMaId(grandmaAmp));
  theMutex.unlock();
  result=_currentParamMap.at(absLamX)*_decAmpDaughter1->XdecAmp(lamX, theData, this)*currentDyn;
  //  result=_currentParamMap.at(absLamX)*_decAmpDaughter1->XdecAmp(lamX, theData, this);
  if(result.real()!=result.real()){
    InfoMsg << "dyn name: " << _absDyn->name() 
	 << "\nname(): " << name()
	 << endmsg;
    Alert << "result:\t" << result << endmsg;
    exit(0);
  }

  if ( _cacheAmps){
    theMutex.lock();
    _cachedAmpIdMap[evtNo][_absDyn->grandMaId(grandmaAmp)][currentSpinIndex]=result;
    theMutex.unlock();
  }
  return result;
}

void FormXDecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return;
  theMutex.lock();
  _cachedDynIdMap[std::this_thread::get_id()][_absDyn->grandMaId(grandmaAmp)] = _absDyn->eval( theData, grandmaAmp);
  theMutex.unlock();
 
  if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 return;
}

void  FormXDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::map<Spin, std::string>::iterator it;
  for (it=_paramNameMap.begin(); it!=_paramNameMap.end(); ++it){
    Spin currentHeli=it->first;
    if( currentHeli >= 0){
      std::string currentName=it->second+"Mag";
      fitPar->Add(currentName, 1., 0.2);
      currentName=it->second+"Phi";
      fitPar->Add(currentName, 0., 0.2);
      _currentParamMap[it->first]=complex<double>(1.,0.);
    }
  }
  _absDyn->fillDefaultParams(fitPar); 
  if(!_daughter1IsStable) _decAmpDaughter1->fillDefaultParams(fitPar);
}

void FormXDecAmps::fillParamNameList(){
  _paramNameList.clear();
  Spin lamRes(0);
  _paramNameMap[lamRes]=absDec()->name()+"lam0";
  _paramNameList.push_back(_paramNameMap.at(lamRes)+"Mag");
  _paramNameList.push_back(_paramNameMap.at(lamRes)+"Phi");
  if(_J> 1){
    lamRes=-2;
    _paramNameMap[lamRes]=absDec()->name()+"lam2";
    _paramNameList.push_back(_paramNameMap.at(lamRes)+"Mag");
    _paramNameList.push_back(_paramNameMap.at(lamRes)+"Phi");
    lamRes=2;
    _paramNameMap[lamRes]=absDec()->name()+"lam2";
    _paramNameList.push_back(_paramNameMap.at(lamRes)+"Mag");
    _paramNameList.push_back(_paramNameMap.at(lamRes)+"Phi");
  }
}

void FormXDecAmps::print(std::ostream& os) const{
  return; //dummy
}


void FormXDecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::map<Spin, std::string>::iterator it;
  for (it=_paramNameMap.begin(); it!=_paramNameMap.end(); ++it){
    Spin currentHeli=it->first;
    if( currentHeli >= 0){
      std::string currentName=it->second+"Mag";
      double theMag=fabs(fitPar->Value(currentName));
      currentName=it->second+"Phi";
      double thePhi=fitPar->Value(currentName);
      complex<double> expi(cos(thePhi), sin(thePhi)); 
      _currentParamMap.at(it->first)=theMag*expi;
    }
  }
  _absDyn->updateFitParams(fitPar);  
  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
}

void FormXDecAmps::initialize(){
  if( _daughter1IsStable){
    Alert << "for the formation amplitude the daughter1 must be an unstable paticle!!!" << endmsg;
    exit(1);
  }

  if(!_daughter1IsStable){
    std::shared_ptr<AbsDecay> decDaughter1=_decay->decDaughter1();
    std::shared_ptr<AbsDynamics> currentDyn=decDaughter1->getDynamics();
    if (0!=currentDyn) decDaughter1->getDynamics()->addGrandMa(_decay);
    else{
      WarningMsg <<"no dynamics for\t" << decDaughter1->name() << endmsg;
    }
    _decAmpDaughter1=XdecAmpRegistry::instance()->getXdecAmp(_channelID, decDaughter1);
  }

  _Jdaughter1=(Spin) _decay->daughter1Part()->J();
   
  _lam1Min = -_Jdaughter1;
  _lam1Max = _Jdaughter1;

  _daughter1Name=_decay->daughter1Part()->name();

}


complex<double> FormXDecAmps::daughterAmp(Spin& lam1, Spin& lam2, EvtData* theData){
  complex<double> result(1.,0.);
  result *= _decAmpDaughter1->XdecAmp(lam1, theData, this);
  return result;
}

void FormXDecAmps::cacheAmplitudes(){
  _cacheAmps=true;
  _decAmpDaughter1->cacheAmplitudes();
}

bool FormXDecAmps::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){
  _recalculate=false;
  if(_decAmpDaughter1->checkRecalculation(fitParNew, fitParOld)) _recalculate=true;
  if (!_recalculate) _recalculate=AbsParamHandler::checkRecalculation(fitParNew, fitParOld);
  if (!_recalculate) _recalculate=_absDyn->checkRecalculation(fitParNew, fitParOld);

  return _recalculate;
}

