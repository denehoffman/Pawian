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

// LSOmegaTo3PiDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/LSOmegaTo3PiDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/OmegaTo3PiLSDecay.hh"
//#include "PwaUtils/XdecAmpRegistry.hh"
#include "Particle/Particle.hh"


LSOmegaTo3PiDecAmps::LSOmegaTo3PiDecAmps(boost::shared_ptr<OmegaTo3PiLSDecay> theDec) :
  LSDecAmps(theDec)
  ,_lambdaDecKey(theDec->lambdaDecKey())
{
  _JPCLSs=theDec->JPCLSAmps();
  _factorMag=1.;
  if(_JPCLSs.size()>0) _factorMag=1./sqrt(_JPCLSs.size());
}

LSOmegaTo3PiDecAmps::~LSOmegaTo3PiDecAmps()
{
}


complex<double> LSOmegaTo3PiDecAmps::XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));

        complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
	  *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][0]);
        result+=amp;
  }

  return result;
}




complex<double> LSOmegaTo3PiDecAmps::XdecAmp(Spin lamX, EvtData* theData, Spin lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.); 
  
  int evtNo=theData->evtNo;
  std::string currentKey=_absDyn->grandMaKey(grandmaAmp); 

  if ( _cacheAmps && !_recalculate){
    result=_cachedGrandmaAmpMap[currentKey][evtNo][lamX][lamFs];
    return result;
  }
  
  std::vector< boost::shared_ptr<const JPCLS> >::iterator it;
  for (it=_JPCLSs.begin(); it!=_JPCLSs.end(); ++it){
    if( fabs(lamX) > (*it)->J ) continue;
    double theMag=_currentParamMags[*it];
    double thePhi=_currentParamPhis[*it];
    complex<double> expi(cos(thePhi), sin(thePhi));
    
    complex<double> amp = theMag*expi*sqrt(2*(*it)->L+1)
      *conj( theData->WignerDsString[_wignerDKey][(*it)->J][lamX][0]);
    
    result+=amp;
  }
  result*=sqrt( theData->DoubleString["lamOmegaDec"] );
  result*=_absDyn->eval(theData, grandmaAmp);

  if ( _cacheAmps){
     theMutex.lock();
     _cachedGrandmaAmpMap[currentKey][evtNo][lamX][lamFs]=result;
     theMutex.unlock();
  }
  return result;
}


void LSOmegaTo3PiDecAmps::print(std::ostream& os) const{
  return; //dummy
}


