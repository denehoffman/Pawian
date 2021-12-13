//************************************************************************//
//									  //
//  Copyright 2021 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// VoigtBlattWRelDynamics class definition file. -*- C++ -*-
// Copyright 20213Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/VoigtBlattWRelDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/BreitWignerFunction.hh"
#include "Utils/IdStringMapRegistry.hh"
#include "Utils/MathUtils.hh"
#include "Utils/Faddeeva.hh"
#include "qft++Extension/PawianUtils.hh"

VoigtBlattWRelDynamics::VoigtBlattWRelDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother, std::vector<Particle*>& fsParticlesDaughter1, std::vector<Particle*>& fsParticlesDaughter2, double qR) :
  BreitWignerBlattWRelDynamics(name, fsParticles, mother, fsParticlesDaughter1, fsParticlesDaughter2, qR)
  ,_massSigmaKey("defaultMassSigmaKey")
  ,_currentSigma(1.)
{
  _isLdependent=true;
}

VoigtBlattWRelDynamics::~VoigtBlattWRelDynamics()
{
}

complex<double> VoigtBlattWRelDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

  int orbMom(OrbMom);

  if (!_recalculate){
    return _cachedLMap.at(theData->evtNo).at(orbMom);
  }
  //the convolution with a Gassian must be put in here

//    double massA = theData->DoubleMassId.at(_dynMassIdDaughter1);
//    double massB = theData->DoubleMassId.at(_dynMassIdDaughter2);
//    
//    complex<double> i(0.,1.);
//    complex<double> rho0=PawianQFT::phaseSpaceFacDefault(_currentMass, massA, massB);
//    complex<double> rho=PawianQFT::phaseSpaceFacDefault(theData->DoubleMassId.at(_dynId), massA, massB);
//    complex<double> momQ0=PawianQFT::breakupMomQDefault(_currentMass, massA, massB);
//    complex<double> momQ=PawianQFT::breakupMomQDefault(theData->DoubleMassId.at(_dynId), massA, massB);
//
//    double width = _currentWidth*std::norm(BarrierFactor::BlattWeisskopfRatio(orbMom, momQ, momQ0, _qR)); 
//
//    double temp=0.;
//    if ((_currentSigma < 0. || width < 0.) || (fabs(_currentSigma)<1e-20 && fabs(width)<1.e-20)) {
//        return temp;  // Not meant to be for those who want to be thinner than 0
//    }
//    if (fabs(_currentSigma) < 1.e-20){
//        _currentSigma=1.e-12;
//    }
//    if (fabs(width) < 1.e-20){
//        width=1.e-12;
//    }
//
//    double denom=sqrt(2.)*_currentSigma;  
//    double realZ=(theData->DoubleMassId.at(_dynId)-_currentMass)/denom;
//    double imagZ=width/(2.*denom);
//    complex<double> complZ(realZ,imagZ);  
//
//    temp=sqrt(2.*M_PI)/4.*_currentWidth*std::abs(BarrierFactor::BlattWeisskopfRatio(orbMom, momQ, momQ0, _qR))/_currentSigma*Faddeeva::w(complZ).real();
//
//    complex<double> result(sqrt(temp), 0.);

    complex<double> result(0.,0.);
    //int nConv = 100;

    double xMin = theData->DoubleMassId.at(_dynId)-3.*_currentWidth;
    double xMax = theData->DoubleMassId.at(_dynId)+3.*_currentWidth;
    
    int nConv = std::ceil(6.*_currentWidth/0.002);
    
    double step = (xMax-xMin)/((double)nConv);
    double mean = 0.; 
    
    
   for(int i=0; i<nConv; i++){
        double xx = xMin+i*step;

        double y1 = std::abs(BreitWignerFunction::BlattWRel(orbMom, xx, _currentMass, _currentWidth, theData->DoubleMassId.at(_dynMassIdDaughter1), theData->DoubleMassId.at(_dynMassIdDaughter2), _qR)) * MathUtils::Gauss(theData->DoubleMassId.at(_dynId)-xx, mean, _currentSigma);

        xx = xMin+(i+1)*step;

        double y2 = std::abs(BreitWignerFunction::BlattWRel(orbMom, xx, _currentMass, _currentWidth, theData->DoubleMassId.at(_dynMassIdDaughter1), theData->DoubleMassId.at(_dynMassIdDaughter2), _qR)) * MathUtils::Gauss(theData->DoubleMassId.at(_dynId)-xx, mean, _currentSigma);

        result+= 0.5*(y2+y1)*step;
    }

  if ( _cacheAmps){
     theMutex.lock();
     _cachedLMap[theData->evtNo][orbMom]=result;
     theMutex.unlock();
  }  

  return result;
}

void  VoigtBlattWRelDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  BreitWignerBlattWRelDynamics::fillDefaultParams(fitPar);
    //fill sigma width
  fitPar->Add(_massSigmaKey, 0.01, 0.4*0.01);
  fitPar->SetLimits(_massSigmaKey, 0., 0.06);
}

void VoigtBlattWRelDynamics::fillParamNameList(){
  _paramNameList.clear();
  BreitWignerBlattWRelDynamics::fillParamNameList();

  //fill sigma width
  _paramNameList.push_back(_massSigmaKey);
}

void VoigtBlattWRelDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  BreitWignerBlattWRelDynamics::updateFitParams(fitPar);
  _currentSigma=fitPar->Value(_massSigmaKey);
}

void VoigtBlattWRelDynamics::setMassKey(std::string& theMassKey){
  BreitWignerBlattWRelDynamics::setMassKey(theMassKey);
  _massSigmaKey=theMassKey+"Sigma";
}
