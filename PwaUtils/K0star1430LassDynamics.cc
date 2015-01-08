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

void  K0star1430LassDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
    fitVal.Masses[_massKey]=_mother->mass();
    fitErr.Masses[_massKey]=0.03;
    fitVal.Widths[_massKey]=_mother->width();
    fitErr.Widths[_massKey]=0.2*_mother->width();
    fitVal.otherParams[_aLASSKey]=1.07;
    fitErr.otherParams[_aLASSKey]=0.05;    
    fitVal.otherParams[_rLASSKey]=-2.852;
    fitErr.otherParams[_rLASSKey]=0.05;  
    fitVal.otherParams[_BLASSKey]=0.7;
    fitErr.otherParams[_BLASSKey]=0.05;  
    fitVal.otherParams[_phiBKey]=0.7;
    fitErr.otherParams[_phiBKey]=0.05;
    fitVal.otherParams[_RLASSKey]=1.;
    fitErr.otherParams[_RLASSKey]=0.05;
    fitVal.otherParams[_phiRKey]=-5.356;
    fitErr.otherParams[_phiRKey]=0.05;
}

bool K0star1430LassDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  double mass=theParamVal.Masses[_massKey];
  if (!CheckDoubleEquality(mass, _currentMass)){
    _recalculate=true;
    return _recalculate;
  }

  double width=theParamVal.Widths[_massKey];
  if (!CheckDoubleEquality(width, _currentWidth)){
    _recalculate=true;
    return _recalculate;
  }

  double aLASS=theParamVal.otherParams[_aLASSKey];
  if (!CheckDoubleEquality(aLASS, _currentaLASS)){
    _recalculate=true;
    return _recalculate;
  }

  double rLASS=theParamVal.otherParams[_rLASSKey];
  if (!CheckDoubleEquality(rLASS, _currentrLASS)){
    _recalculate=true;
    return _recalculate;
  }

  double BLASS=theParamVal.otherParams[_BLASSKey];
  if (!CheckDoubleEquality(BLASS, _currentBLASS)){
    _recalculate=true;
    return _recalculate;
  }

  double phiB=theParamVal.otherParams[_phiBKey];
  if (!CheckDoubleEquality(phiB, _currentphiB)){
    _recalculate=true;
    return _recalculate;
  }

  double RLASS=theParamVal.otherParams[_RLASSKey];
  if (!CheckDoubleEquality(RLASS, _currentRLASS)){
    _recalculate=true;
    return _recalculate;
  }

  double phiR=theParamVal.otherParams[_phiRKey];
  if (!CheckDoubleEquality(phiR, _currentphiR)){
    _recalculate=true;
    return _recalculate;
  }

  return _recalculate;
}

void K0star1430LassDynamics::updateFitParams(fitParams& theParamVal){
  _currentMass=theParamVal.Masses.at(_massKey);
  _currentWidth=theParamVal.Widths.at(_massKey);
  _currentaLASS=theParamVal.otherParams.at(_aLASSKey);
  _currentrLASS=theParamVal.otherParams.at(_rLASSKey);
  _currentBLASS=theParamVal.otherParams.at(_BLASSKey);
  _currentphiB=theParamVal.otherParams.at(_phiBKey);
  _currentRLASS=theParamVal.otherParams.at(_RLASSKey);
  _currentphiR=theParamVal.otherParams.at(_phiRKey);
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
