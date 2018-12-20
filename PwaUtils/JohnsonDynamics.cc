//************************************************************************//
//                                                                        //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)                      //
//               Julian Pychy (julian@ep1.rub.de)                         //
//               - Ruhr-Universität Bochum                                //
//                                                                        //
//  This file is part of Pawian.                                          //
//                                                                        //
//  Pawian is free software: you can redistribute it and/or modify        //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or     //
//  (at your option) any later version.                                   //
//                                                                        //
//  Pawian is distributed in the hope that it will be useful,             //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of        //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
//  GNU General Public License for more details.                          //
//                                                                        //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.       //
//                                                                        //
//************************************************************************//

// JohnsonDynamics class definition file. -*- C++ -*-
// Copyright 2018 Andreas Pitka

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>
using namespace std;

#include "PwaUtils/JohnsonDynamics.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/Johnson.hh"
#include "FitParams/AbsPawianParameters.hh"

JohnsonDynamics::JohnsonDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother),
  _currentMass0(1.),
  _currentSigma(1.),
  _currentDelta(1.),
  _currentGamma(1.),
  _johnsonPtr(new Johnson())
{
  _isLdependent=false;
}

JohnsonDynamics::~JohnsonDynamics()
{
}

complex<double> JohnsonDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){
  //  int evtNo=theData->evtNo;
  if (!_recalculate){
    return _cachedMap[theData->evtNo];
  }

  complex<double> result( sqrt( _johnsonPtr->calc(theData->DoubleMassId.at(_dynId), _currentMass0, _currentSigma, _currentDelta, _currentGamma) ), 0.);
  
  if ( _cacheAmps){
     theMutex.lock();
     _cachedMap[theData->evtNo]=result;
     theMutex.unlock();
  }  

  return result;
}

void  JohnsonDynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  // fill mass
  std::string massName=_massKey+"Mass";
  double valMass=_mother->mass();
  double errMass=0.03;
  double minMass=valMass-5.*errMass;
  if(minMass<0.) minMass=0.;
  double maxMass=valMass+5.*errMass;
  
  fitPar->Add(massName, valMass, errMass);
  fitPar->SetLimits(massName, minMass, maxMass);

  // fill sigma   
  std::string sigmaName=_massKey+"Sigma";
  fitPar->Add(sigmaName, 0.44e-3, 0.44e-3*0.01);
  fitPar->SetLimits(sigmaName, 0., 0.01);

  // fill delta   
  std::string deltaName=_massKey+"Delta";
  fitPar->Add(deltaName, 0.69, 0.69*0.01);
  fitPar->SetLimits(deltaName, -5., 5.);

  // fill gamma 
  std::string gammaName=_massKey+"Gamma";
  fitPar->Add(gammaName, -0.15, 0.15*0.01);
  fitPar->SetLimits(gammaName, -10., 10.);

}

void JohnsonDynamics::fillParamNameList(){
  _paramNameList.clear();

  // mass
  std::string massName=_massKey+"Mass";
  _paramNameList.push_back(massName);

  // sigma
  std::string sigmaName=_massKey+"Sigma";
  _paramNameList.push_back(sigmaName);

  // delta
  std::string deltaName=_massKey+"Delta";
  _paramNameList.push_back(deltaName);

  // gamma
  std::string gammaName=_massKey+"Gamma";
  _paramNameList.push_back(gammaName);

}

void JohnsonDynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::string massName=_massKey+"Mass";
  _currentMass0=fitPar->Value(massName);

  std::string sigmaName=_massKey+"Sigma";
  _currentSigma=fitPar->Value(sigmaName);

  std::string deltaName=_massKey+"Delta";
  _currentDelta=fitPar->Value(deltaName);

  std::string gammaName=_massKey+"Gamma";
  _currentGamma=fitPar->Value(gammaName);
}

void JohnsonDynamics::setMassKey(std::string& theMassKey){
  _massKey=theMassKey;
}
