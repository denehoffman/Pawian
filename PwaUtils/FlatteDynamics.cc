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

// FlatteDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>
#include <boost/algorithm/string.hpp>

#include "PwaUtils/FlatteDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/Flatte.hh"
#include "Utils/FunctionUtils.hh"

FlatteDynamics::FlatteDynamics(std::string& key, std::vector<Particle*>& fsParticles, Particle* mother, std::pair<Particle*, Particle*>& decPair1stChannel, std::pair<Particle*, Particle*>& decPair2ndChannel) :
  AbsDynamics(key, fsParticles, mother)
  ,_g11Key(_massKey+decPair1stChannel.first->name()+decPair1stChannel.second->name())
  ,_g22Key(_massKey+decPair2ndChannel.first->name()+decPair2ndChannel.second->name())
{
  std::string firstPair1stPartStr=decPair1stChannel.first->name();
  boost::replace_all(firstPair1stPartStr,"+", "");
  boost::replace_all(firstPair1stPartStr,"-", "");

  std::string firstPair2ndPartStr=decPair1stChannel.second->name();
  boost::replace_all(firstPair2ndPartStr,"+", "");
  boost::replace_all(firstPair2ndPartStr,"-", "");

  _g11Key=_massKey+firstPair1stPartStr+firstPair2ndPartStr;

  std::string secondPair1stPartStr=decPair2ndChannel.first->name();
  boost::replace_all(secondPair1stPartStr,"+", "");
  boost::replace_all(secondPair1stPartStr,"-", "");

  std::string secondPair2ndPartStr=decPair2ndChannel.second->name();
  boost::replace_all(secondPair2ndPartStr,"+", "");
  boost::replace_all(secondPair2ndPartStr,"-", "");

  _g22Key=_massKey+secondPair1stPartStr+secondPair2ndPartStr;  

  _flattePtr=std::shared_ptr<Flatte>(new Flatte(decPair1stChannel, decPair2ndChannel));
}

FlatteDynamics::~FlatteDynamics()
{
}

complex<double> FlatteDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result=_flattePtr->calcFirstChannel(theData->FourVecsString[_dynKey].M(), _currentMass, _currentg11, _currentg22);
  
  if ( _cacheAmps){
      theMutex.lock();
      _cachedMap[evtNo]=result;
      theMutex.unlock();
  }  

  return result;
}

void  FlatteDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
    fitVal.Masses[_massKey]=_mother->mass();
    fitErr.Masses[_massKey]=0.03;
    fitVal.gFactors[_g11Key]=1.;
    fitErr.gFactors[_g11Key]=1.;
    fitVal.gFactors[_g22Key]=1.;
    fitErr.gFactors[_g22Key]=1.;
}

bool FlatteDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  double mass=theParamVal.Masses[_massKey];
  if (!CheckDoubleEquality(mass, _currentMass)){
    _recalculate=true;
  }

  double g11=theParamVal.gFactors[_g11Key];
  if (!CheckDoubleEquality(g11, _currentg11)){
    _recalculate=true;
  }

  double g22=theParamVal.gFactors[_g22Key];
  if (!CheckDoubleEquality(g22, _currentg22)){
    _recalculate=true;
  }
  return _recalculate;
}

void FlatteDynamics::updateFitParams(fitParams& theParamVal){
  _currentMass=theParamVal.Masses[_massKey];
  _currentg11=theParamVal.gFactors[_g11Key];
  _currentg22=theParamVal.gFactors[_g22Key];
}

void FlatteDynamics::setMassKey(std::string& theMassKey){
  boost::replace_all(_g11Key,_massKey, theMassKey);
  boost::replace_all(_g22Key,_massKey, theMassKey);
  AbsDynamics::setMassKey(theMassKey);
}
