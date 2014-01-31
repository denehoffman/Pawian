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

// TensorDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <mutex>

#include "PwaUtils/TensorDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/IsobarTensorDecay.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"

TensorDecAmps::TensorDecAmps(std::shared_ptr<IsobarTensorDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_LSs(theDec->LSAmps())
  ,_factorMag(1.)
{
  if(_LSs.size()>0) _factorMag=1./sqrt(_LSs.size());
}

TensorDecAmps::TensorDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
{
  if(_LSs.size()>0) _factorMag=1./sqrt(_LSs.size());
}

TensorDecAmps::~TensorDecAmps()
{
}


complex<double> TensorDecAmps::XdecPartAmp(Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  Spin lam1Min=-_Jdaughter1;
  Spin lam1Max= _Jdaughter1;
  Spin lam2Min=-_Jdaughter2;
  Spin lam2Max=_Jdaughter2;

  if(fixDaughterNr == 1){
     lam1Min = lam1Max = lamDec;
  }
  else if(fixDaughterNr == 2){
     lam2Min = lam2Max = lamDec;
  }
  else{
     Alert << "Invalid fixDaughterNr in XdecPartAmp." << endmsg;
  }

  if(_enabledlamFsDaughter1){
    lam1Min=lamFs;
    lam1Max=lamFs;
  }
  else if(_enabledlamFsDaughter2){
    lam2Min=lamFs;
    lam2Max=lamFs;
  }

  complex<double> result=lsLoop(lamX, theData, lam1Min, lam1Max, lam2Min, lam2Max, false);

  return result;
}




complex<double> TensorDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;

  if ( _cacheAmps && !_recalculate){
    result=_cachedAmpMap.at(evtNo).at(lamX).at(lamFs);
    result*=_absDyn->eval(theData, grandmaAmp);
    return result;
  }

  //  Spin lam1Min=-_Jdaughter1;
  Spin lam1Min=-_Jdaughter1;
  Spin lam1Max= _Jdaughter1;
  Spin lam2Min=-_Jdaughter2;
  Spin lam2Max=_Jdaughter2;

  if(_enabledlamFsDaughter1){
    lam1Min=lamFs;
    lam1Max=lamFs;
  }
  else if(_enabledlamFsDaughter2){
    lam2Min=lamFs;
    lam2Max=lamFs;
  }


  result=lsLoop(lamX, theData, lam1Min, lam1Max, lam2Min, lam2Max, true, lamFs );

  if ( _cacheAmps){
     theMutex.lock();
     _cachedAmpMap[evtNo][lamX][lamFs]=result;
     theMutex.unlock();
  }
  // if (norm(result)>1.){
  //   Info << "evtNo: " << evtNo << "  name(): " << name() << " result(lamX=" << lamX << ", lamFs=" << lamFs << "): " << result << endmsg;
  // }
  result*=_absDyn->eval(theData, grandmaAmp);
  return result;
}


complex<double> TensorDecAmps::lsLoop(Spin lamX, EvtData* theData, Spin lam1Min, Spin lam1Max, Spin lam2Min, Spin lam2Max, bool withDecs, Spin lamFs ){
  complex<double> result(0.,0.);
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    double theMag=_currentParamMags.at(*it);
    double thePhi=_currentParamPhis.at(*it);
    complex<double> expi(cos(thePhi), sin(thePhi));

    for(Spin lambda1=lam1Min; lambda1<=lam1Max; ++lambda1){
      for(Spin lambda2=lam2Min; lambda2<=lam2Max; ++lambda2){
	//	Spin lambda = lambda1-lambda2;
	//	if((*it)->S < fabs(lambda) || _JPCPtr->J < fabs(lambda)) continue;
       complex<double> amp = theMag*expi*theData->ComplexDouble5SpinString.at(_name).at((*it)->L).at((*it)->S).at(lamX).at(lambda1).at(lambda2);
       // if (norm(amp)>1000.){
       //   Info << "norm(amp): " << norm(amp) << endmsg;
       // 	 Info << "evtNo: " << theData->evtNo << "  name(): " << name() << " result(L=" << (*it)->L << ", S=" << (*it)->S << ", lamX=" << lamX << ", lambda1=" << lambda1 << ", lambda2=" << lambda2 << "): " << amp << endmsg;
       // 	 //	 amp= complex<double>(0.,0.);
       // }
      	if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData, lamFs);
	result+=amp;
      }
    }
  }

  result*=_isospinCG;
  return result;
}


void  TensorDecAmps::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
  std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    currentMagValMap[*itLS]=_factorMag;
    currentPhiValMap[*itLS]=0.;
    currentMagErrMap[*itLS]=_factorMag;
    currentPhiErrMap[*itLS]=0.3;
  }

  fitVal.MagsLS[_key]=currentMagValMap;
  fitVal.PhisLS[_key]=currentPhiValMap;
  fitErr.MagsLS[_key]=currentMagErrMap;
  fitErr.PhisLS[_key]=currentPhiErrMap;

  _absDyn->getDefaultParams(fitVal, fitErr);


  if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
  if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);
}

void TensorDecAmps::print(std::ostream& os) const{
  return; //dummy
}


bool TensorDecAmps::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

   if(_absDyn->checkRecalculation(theParamVal)) _recalculate=true;

   if(!_daughter1IsStable) {
     if(_decAmpDaughter1->checkRecalculation(theParamVal)) _recalculate=true;
   }
   if(!_daughter2IsStable){
     if(_decAmpDaughter2->checkRecalculation(theParamVal)) _recalculate=true;
   }

   if(!_recalculate){
     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagsLS[_key];
     std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhisLS[_key];

     std::vector< std::shared_ptr<const LScomb> >::iterator it;
     for (it=_LSs.begin(); it!=_LSs.end(); ++it){
       double theMag=magMap[*it];
       double thePhi=phiMap[*it];

       if(!CheckDoubleEquality(theMag, _currentParamMags[*it])){
	 _recalculate=true;
	 return _recalculate;
       }
       if(!CheckDoubleEquality(thePhi, _currentParamPhis[*it])){
	 _recalculate=true;
	 return _recalculate;
       }
     }
   }

   return _recalculate;
}


void  TensorDecAmps::updateFitParams(fitParams& theParamVal){
   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& magMap=theParamVal.MagsLS[_key];
   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& phiMap=theParamVal.PhisLS[_key];

   std::vector< std::shared_ptr<const LScomb> >::iterator it;
   for (it=_LSs.begin(); it!=_LSs.end(); ++it){
     double theMag=magMap[*it];
     double thePhi=phiMap[*it];
     _currentParamMags[*it]=theMag;
     _currentParamPhis[*it]=thePhi;
   }

   _absDyn->updateFitParams(theParamVal);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(theParamVal);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(theParamVal);

}

