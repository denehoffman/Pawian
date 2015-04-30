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
#include "FitParams/AbsPawianParameters.hh"

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
  _isLdependent=false;
}

FlatteDynamics::~FlatteDynamics()
{
}

complex<double> FlatteDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result=_flattePtr->calcFirstChannel(theData->DoubleString.at(_dynKey), _currentMass, _currentg11, _currentg22);
  
  if ( _cacheAmps){
      theMutex.lock();
      _cachedMap[evtNo]=result;
      theMutex.unlock();
  }  

  return result;
}

void FlatteDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill mass
  std::string massName=_massKey+"Mass";
  double valMass=_mother->mass();
  double errMass=0.03;
  double minMass=valMass-5.*errMass;
  if(minMass<0.) minMass=0.;
  double maxMass=valMass+5.*errMass;
  
  fitPar->Add(massName, valMass, errMass);
  fitPar->SetLimits(massName, minMass, maxMass);
  
  //fill g-Factors
  std::string g11Name=_g11Key+"gFactor";
  double valgii=1.;
  double errgii=0.5;
  double mingii=valgii-5.*errgii;
  double maxgii=valgii+20.*errgii;
  fitPar->Add(g11Name, valgii, errgii);
  fitPar->SetLimits(g11Name, mingii, maxgii);

  std::string g22Name=_g22Key+"gFactor";
  fitPar->Add(g22Name, valgii, errgii);
  fitPar->SetLimits(g22Name, mingii, maxgii);
}

void FlatteDynamics::fillParamNameList(){
 std::string massName=_massKey+"Mass";
  _paramNameList.push_back(massName);  
  
  //g-Factors
  std::string g11Name=_g11Key+"gFactor";
  _paramNameList.push_back(g11Name);  

  std::string g22Name=_g22Key+"gFactor";
  _paramNameList.push_back(g22Name);  
}

void  FlatteDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::string massName=_massKey+"Mass";
  _currentMass=fitPar->Value(massName);
  
  //g-Factors
  std::string g11Name=_g11Key+"gFactor";
  _currentg11=fitPar->Value(g11Name);

  std::string g22Name=_g22Key+"gFactor";
  _currentg22=fitPar->Value(g22Name);  
}

void FlatteDynamics::setMassKey(std::string& theMassKey){
  boost::replace_all(_g11Key,_massKey, theMassKey);
  boost::replace_all(_g22Key,_massKey, theMassKey);
  AbsDynamics::setMassKey(theMassKey);
}
