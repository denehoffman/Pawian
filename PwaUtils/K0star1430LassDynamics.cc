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

// K0star1430LassDynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>
#include <boost/algorithm/string.hpp>

#include "PwaUtils/K0star1430LassDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/LASS.hh"
#include "Utils/FunctionUtils.hh"
#include "FitParams/AbsPawianParameters.hh"

K0star1430LassDynamics::K0star1430LassDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
  ,_aLASSKey(_massKey+"aLASS_PosNeg")
  ,_rLASSKey(_massKey+"rLASS_PosNeg")
  ,_BLASSKey(_massKey+"BLASS_PosNeg")
  ,_phiBKey(_massKey+"phiB_PosNeg")
  ,_RLASSKey(_massKey+"RLASS_PosNeg")
  ,_phiRKey(_massKey+"phiR_PosNeg")
{
  _isLdependent=false;
}

K0star1430LassDynamics::~K0star1430LassDynamics()
{
}

complex<double> K0star1430LassDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result=LASS::K0star_1430(theData->DoubleString.at(_dynKey), _currentMass, _currentWidth, _currentaLASS, _currentrLASS, _currentBLASS, _currentphiB, _currentRLASS, _currentphiR);

  
  if ( _cacheAmps){
      theMutex.lock();
      _cachedMap[evtNo]=result;
      theMutex.unlock();
  }  

  return result;
}

void  K0star1430LassDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string massName=_massKey+"Mass";
  double valMass=_mother->mass();
  double errMass=0.03;
  double minMass=valMass-5.*errMass;
  if(minMass<0.) minMass=0.;
  double maxMass=valMass+5.*errMass;
  fitPar->Add( massName, valMass, errMass);
  fitPar->SetLimits( massName, minMass, maxMass);

  //fill width
  std::string widthName=_massKey+"Width";
  double valWidth=_mother->width();
  double errWidth=0.2*_mother->width();
  double minWidth=valWidth-5.*errWidth;
  if(minWidth<0.) minWidth=0.;
  double maxWidth=valWidth+5.*errWidth;

  fitPar->Add(widthName, valWidth, errWidth);
  fitPar->SetLimits(widthName, minWidth, maxWidth);

  fitPar->Add(_aLASSKey, 1.07, 0.05);
  fitPar->Add(_rLASSKey, -2.852, 0.05);
  fitPar->Add(_BLASSKey, 0.7, 0.05);
  fitPar->Add(_phiBKey, 0.7, 0.05);
  fitPar->Add(_RLASSKey, 1., 0.05);
  fitPar->Add(_phiRKey, -5.356, 0.05);
}

void K0star1430LassDynamics::fillParamNameList(){
  _paramNameList.clear();  

  std::string massName=_massKey+"Mass";
  _paramNameList.push_back(massName);

  std::string widthName=_massKey+"Width";
  _paramNameList.push_back(widthName);

  _paramNameList.push_back(_aLASSKey);
  _paramNameList.push_back(_rLASSKey);
  _paramNameList.push_back(_BLASSKey);
  _paramNameList.push_back(_RLASSKey);
  _paramNameList.push_back(_phiRKey);
} 

void K0star1430LassDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string massName=_massKey+"Mass";
  _currentMass=fitPar->Value(massName);

  std::string widthName=_massKey+"Width";
  _currentWidth=fitPar->Value(widthName);

  _currentaLASS=fitPar->Value(_aLASSKey);
  _currentrLASS=fitPar->Value(_rLASSKey);
  _currentBLASS=fitPar->Value(_BLASSKey);
  _currentphiB=fitPar->Value(_phiBKey);
  _currentRLASS=fitPar->Value(_RLASSKey);
  _currentphiR=fitPar->Value(_phiRKey);

}

void K0star1430LassDynamics::setMassKey(std::string& theMassKey){
  boost::replace_all(_aLASSKey,_massKey, theMassKey);
  boost::replace_all(_rLASSKey,_massKey, theMassKey);
  boost::replace_all(_BLASSKey,_massKey, theMassKey);
  boost::replace_all(_phiBKey,_massKey, theMassKey);
  boost::replace_all(_RLASSKey,_massKey, theMassKey);
  boost::replace_all(_phiRKey,_massKey, theMassKey);
  AbsDynamics::setMassKey(theMassKey);
}
