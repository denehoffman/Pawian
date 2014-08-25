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

// AbsXdecAmp class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/DynRegistry.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"


AbsXdecAmp::AbsXdecAmp(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsParamHandler()
  , _channelID(channelID)
  , _decay(theDec)
  , _name(theDec->name())
  ,_JPCPtr(theDec->motherJPC())
  , _jpcDecsName(_JPCPtr->jpcname()+"To"+theDec->daughter1Part()->name()+"_"+theDec->daughter2Part()->name())
  , _isoKey(_JPCPtr->jpcname()+theDec->fitParSuffix())
  ,_absDyn(theDec->getDynamics())
  ,_isospinCG(theDec->isospinCG())
  ,_preFactor(theDec->preFactor())
  ,_key("_"+theDec->fitParSuffix())
  ,_wignerDKey(theDec->wignerDKey())
  ,_refKey(theDec->refKey())
  ,_daughter1IsStable(theDec->isDaughter1Stable())
  ,_daughter2IsStable(theDec->isDaughter2Stable())
  ,_enabledlamFsDaughter1(false)
  ,_enabledlamFsDaughter2(false)
  ,_J(theDec->motherJPC()->J)
{
  initialize();
}

AbsXdecAmp::~AbsXdecAmp()
{
}

void AbsXdecAmp::initialize(){
  if(!_decay->useIsospin()) _isospinCG=1.;

  _absDyn = DynRegistry::instance()->getDynamics(_decay);

  if(!_daughter1IsStable){
    std::shared_ptr<AbsDecay> decDaughter1=_decay->decDaughter1();
    std::shared_ptr<AbsDynamics> currentDyn=decDaughter1->getDynamics();
    if (0!=currentDyn) decDaughter1->getDynamics()->addGrandMa(_decay);
    else{
      Warning <<"no dynamics for\t" << decDaughter1->name() << endmsg;
    }
    _decAmpDaughter1=XdecAmpRegistry::instance()->getXdecAmp(_channelID, decDaughter1);
  }

  if(!_daughter2IsStable){
    std::shared_ptr<AbsDecay> decDaughter2=_decay->decDaughter2();

    std::shared_ptr<AbsDynamics> currentDyn=decDaughter2->getDynamics();
    if (0!=currentDyn) decDaughter2->getDynamics()->addGrandMa(_decay);
    else{
      Warning <<"no dynamics for\t" << decDaughter2->name() << endmsg;
    }
    _decAmpDaughter2=XdecAmpRegistry::instance()->getXdecAmp(_channelID, decDaughter2);
  }

  _Jdaughter1=(Spin) _decay->daughter1Part()->J();
  _Jdaughter2=(Spin) _decay->daughter2Part()->J();

  if( _daughter1IsStable && _Jdaughter1>0) _enabledlamFsDaughter1=true;
  if( _daughter2IsStable && _Jdaughter2>0) _enabledlamFsDaughter2=true;
}

complex<double> AbsXdecAmp::daughterAmp(Spin& lam1, Spin& lam2, EvtData* theData, Spin& lamFs){
  complex<double> result(1.,0.);
  if(!_daughter1IsStable) result *= _decAmpDaughter1->XdecAmp(lam1, theData, lamFs, this);
  if(!_daughter2IsStable) result *= _decAmpDaughter2->XdecAmp(lam2, theData, lamFs, this);
  return result;
}

void AbsXdecAmp::cacheAmplitudes(){
  _cacheAmps=true;
  _absDyn->cacheAmplitudes();
  if(!_daughter1IsStable) _decAmpDaughter1->cacheAmplitudes();
  if(!_daughter2IsStable) _decAmpDaughter2->cacheAmplitudes();
}

