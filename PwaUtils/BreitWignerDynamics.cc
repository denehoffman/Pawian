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

// BreitWignerDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/BreitWignerDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"

BreitWignerDynamics::BreitWignerDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
{
  _isLdependent=false;
}

BreitWignerDynamics::~BreitWignerDynamics()
{
}

complex<double> BreitWignerDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }

  complex<double> result=BreitWignerFunction::NonRel(theData->DoubleString.at(_dynKey), _currentMass, _currentWidth);  
  if ( _cacheAmps){
     theMutex.lock();
     _cachedMap[evtNo]=result;
     theMutex.unlock();
  }  

  return result;
}

void  BreitWignerDynamics::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){
    fitVal.Masses[_massKey]=_mother->mass();
    fitErr.Masses[_massKey]=0.03;
    fitVal.Widths[_massKey]=_mother->width();
    fitErr.Widths[_massKey]=0.2*_mother->width();
}

void BreitWignerDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  //fill mass
  std::string massName=_massKey+"Mass";
  double valMass=_mother->mass();
  double errMass=0.03;
  double minMass=valMass-5.*errMass;
  if(minMass<0.) minMass=0.;
  double maxMass=valMass+5.*errMass;
  
  fitPar->Add(massName, valMass, errMass);
  fitPar->SetLimits(massName, minMass, maxMass);
  
  //fill width
  std::string widthName=_massKey+"Width";
  double valWidth=_mother->width();
  double errWidth=0.2*_mother->width();
  double minWidth=valWidth-5.*errWidth;
  if(minWidth<0.) minWidth=0.;
  double maxWidth=valWidth+5.*errWidth;

  fitPar->Add(widthName, valWidth, errWidth);
  fitPar->SetLimits(widthName, minWidth, maxWidth);
}


bool BreitWignerDynamics::checkRecalculation(fitParCol& theParamVal){
  _recalculate=false;

  double mass=theParamVal.Masses[_massKey];
  if (!CheckDoubleEquality(mass, _currentMass)){
    _recalculate=true;
  }
  double width=theParamVal.Widths[_massKey];
  if (!CheckDoubleEquality(width, _currentWidth)){
    _recalculate=true;
  }

  return _recalculate;
}

void BreitWignerDynamics::updateFitParams(fitParCol& theParamVal){
  _currentMass=theParamVal.Masses[_massKey];
  _currentWidth=theParamVal.Widths[_massKey];
}
