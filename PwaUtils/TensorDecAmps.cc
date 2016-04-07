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
  if(_LSs.size()>0) _factorMag/=sqrt(_LSs.size());
}

TensorDecAmps::TensorDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_factorMag(1.)
{
  if(_LSs.size()>0) _factorMag/=sqrt(_LSs.size());
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

  complex<double> result=lsLoop(grandmaAmp, lamX, theData, lam1Min, lam1Max, lam2Min, lam2Max, false);

  return result;
}




complex<double> TensorDecAmps::XdecAmp(Spin& lamX, EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  int evtNo=theData->evtNo;
  Id2StringType currentSpinIndex=FunctionUtils::spin2Index(lamX,lamFs);
  
  if ( _cacheAmps && !_recalculate){
    //result=_cachedAmpMap.at(evtNo).at(_absDyn->grandMaKey(grandmaAmp)).at(currentSpinIndex);
    result=_cachedAmpIdMap.at(evtNo).at(_absDyn->grandMaId(grandmaAmp)).at(currentSpinIndex);
    return result;
  }

  if(_enabledlamFsDaughter1) result=lsLoop(grandmaAmp, lamX, theData, lamFs, lamFs, _lam2Min, _lam2Max, true, lamFs );
  else if(_enabledlamFsDaughter2) result=lsLoop(grandmaAmp, lamX, theData, _lam1Min, _lam1Max, lamFs, lamFs, true, lamFs ); 
  else result=lsLoop(grandmaAmp, lamX, theData, _lam1Min, _lam1Max, _lam2Min, _lam2Max, true, lamFs );

  if ( _cacheAmps){
     theMutex.lock();
     //     _cachedAmpMap[evtNo][_absDyn->grandMaKey(grandmaAmp)][currentSpinIndex]=result;
     _cachedAmpIdMap[evtNo][_absDyn->grandMaId(grandmaAmp)][currentSpinIndex]=result;
     theMutex.unlock();
  }
  return result;
}


complex<double> TensorDecAmps::lsLoop(AbsXdecAmp* grandmaAmp, Spin lamX, EvtData* theData, Spin lam1Min, Spin lam1Max, Spin lam2Min, Spin lam2Max, bool withDecs, Spin lamFs ){
  complex<double> result(0.,0.);

  map<unsigned short, map<Id3StringType, complex<double> > >& currentLS3SpinMap=theData->ComplexLS3Spin.at(_decay->nameId());
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    map<Id3StringType, complex<double> >& current3SpinMap = currentLS3SpinMap.at((*it)->idnumberLS);
    complex<double> theMagExpi=_currentParamMagExpi.at(*it);

    complex<double> tmpResult(0.,0.);
    for(Spin lambda1=lam1Min; lambda1<=lam1Max; ++lambda1){
      for(Spin lambda2=lam2Min; lambda2<=lam2Max; ++lambda2){
	Id3StringType IdLamXLam1Lam2=FunctionUtils::spin3Index(lamX, lambda1, lambda2);
	complex<double> amp = theMagExpi*current3SpinMap.at(IdLamXLam1Lam2);
      	if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData, lamFs);
	tmpResult+=amp;
      }
    }
    if (_absDyn->isLdependent())  tmpResult*=_cachedDynLSMap.at(std::this_thread::get_id()).at((*it)->L);
    result+=tmpResult;
  }
  
  //if (!_absDyn->isLdependent()) result *=_cachedDynMap.at(std::this_thread::get_id()).at(_absDyn->grandMaKey(grandmaAmp));
  if (!_absDyn->isLdependent()) result *=_cachedDynIdMap.at(std::this_thread::get_id()).at(_absDyn->grandMaId(grandmaAmp));
  
  result*=_isospinCG;
  return result;
}


// void  TensorDecAmps::getDefaultParams(fitParCol& fitVal, fitParCol& fitErr){

//   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
//   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
//   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
//   std::map< std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

//   std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
//   for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
//     currentMagValMap[*itLS]=_factorMag;
//     currentPhiValMap[*itLS]=0.;
//     currentMagErrMap[*itLS]=_factorMag/3.;
//     currentPhiErrMap[*itLS]=0.3;
//   }

//   fitVal.MagsLS[_key]=currentMagValMap;
//   fitVal.PhisLS[_key]=currentPhiValMap;
//   fitErr.MagsLS[_key]=currentMagErrMap;
//   fitErr.PhisLS[_key]=currentPhiErrMap;

//   _absDyn->getDefaultParams(fitVal, fitErr);


//   if(!_daughter1IsStable) _decAmpDaughter1->getDefaultParams(fitVal, fitErr);
//   if(!_daughter2IsStable) _decAmpDaughter2->getDefaultParams(fitVal, fitErr);
// }


void  TensorDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

   std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double valMag=_factorMag;
    double errMag=_factorMag/2.;
    //    double minMag=0.;
    //    double maxMag=_factorMag+30.*errMag;

    fitPar->Add(magName, valMag, errMag);
    //    fitPar->SetLimits(magName, minMag, maxMag);

    std::string phiName=(*itLS)->name()+_key+"Phi";
    double valPhi=0.;
    double errPhi=0.2;
    //no limits for phi parameter
    fitPar->Add(phiName, valPhi, errPhi);
  }

  _absDyn->fillDefaultParams(fitPar);

  if(!_daughter1IsStable) _decAmpDaughter1->fillDefaultParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->fillDefaultParams(fitPar);
}

void TensorDecAmps::fillParamNameList(){
  _paramNameList.clear(); 
  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    std::string magName=(*itLS)->name()+_key+"Mag";
    _paramNameList.push_back(magName);

    std::string phiName=(*itLS)->name()+_key+"Phi";
    _paramNameList.push_back(phiName);
  }
}

void TensorDecAmps::print(std::ostream& os) const{
  return; //dummy
}

void TensorDecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    std::string phiName=(*itLS)->name()+_key+"Phi";
    double theMag=fabs(fitPar->Value(magName));
    double thePhi=fitPar->Value(phiName);

    _currentParamMags[*itLS]=theMag;
    _currentParamPhis[*itLS]=thePhi;
    
    complex<double> expi(cos(thePhi), sin(thePhi));
    _currentParamMagExpi[*itLS]=theMag*expi;
  }

  _absDyn->updateFitParams(fitPar);
  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}


void TensorDecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return; 

  if(!_absDyn->isLdependent()){
    AbsXdecAmp::calcDynamics(theData, grandmaAmp);
    return;
  }

 std::vector< std::shared_ptr<const LScomb> >::iterator it;
 for (it=_LSs.begin(); it!=_LSs.end(); ++it){
   theMutex.lock();
   _cachedDynLSMap[std::this_thread::get_id()][(*it)->L]=_absDyn->eval(theData, grandmaAmp, (*it)->L);
   theMutex.unlock();
 }  

 if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
 return;
}


