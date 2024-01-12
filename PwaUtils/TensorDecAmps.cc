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
  initialize();
  if(_LSs.size()>0) _factorMag/=sqrt(_LSs.size());
}

TensorDecAmps::TensorDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_factorMag(1.)
{
  initialize();
  if(_LSs.size()>0) _factorMag/=sqrt(_LSs.size());
}

TensorDecAmps::~TensorDecAmps()
{
}


complex<double> TensorDecAmps::XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr, EvtData* theData, Spin& lamFs,AbsXdecAmp* grandmaAmp){

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




complex<double> TensorDecAmps::XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( std::abs(lamX) > _JPCPtr->J) return result;

  short currentSpinIndex=FunctionUtils::spin1IdIndex(_projId,lamX);  
  
  if (!_recalculate){
    return _cachedAmpIdMap.at(theData->evtNo).at(_absDyn->grandMaId(grandmaAmp)).at(currentSpinIndex);
  }

  result=lsLoop(grandmaAmp, lamX, theData, _lam1MinProj,
		_lam1MaxProj, 
		_lam2MinProj, 
		_lam2MaxProj, true);

  if ( _cacheAmps){
     _cachedAmpIdMap[theData->evtNo][_absDyn->grandMaId(grandmaAmp)][currentSpinIndex]=result;
  }
  return result;
}


complex<double> TensorDecAmps::lsLoop(AbsXdecAmp* grandmaAmp, Spin lamX, EvtData* theData, Spin lam1Min, 
				      Spin lam1Max, Spin lam2Min, Spin lam2Max, bool withDecs) {
  complex<double> result(0.,0.);

  // map<unsigned short, map<Id3StringType, complex<double> > >& currentLS3SpinMap=theData->ComplexLS3Spin.at(_decay->nameId());
  map<Id3StringType, map<unsigned short, complex<double> > >& current3SpinLSMap=theData->Complex3SpinLS.at(_decay->nameId());
  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;

  for(Spin lambda1=-_Jdaughter1; lambda1<=_Jdaughter1; ++lambda1){
    if(lambda1<lam1Min || lambda1>lam1Max) continue;
      for(Spin lambda2=-_Jdaughter2; lambda2<=_Jdaughter2; ++lambda2){ 
      if(lambda2<lam2Min || lambda2>lam2Max) continue;
      Id3StringType IdLamXLam1Lam2=FunctionUtils::spin3Index(lamX, lambda1, lambda2);
      map<unsigned short, complex<double> >& currentLSMap = current3SpinLSMap.at(IdLamXLam1Lam2);
      complex<double> amp(0.,0.);
      for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
	//complex<double> tmpamp=_currentParamMagExpi.at(*itLS)*currentLS3SpinMap.at((*itLS)->idnumberLS).at(IdLamXLam1Lam2);
	complex<double> tmpamp=_currentParamMagExpi.at(*itLS)*currentLSMap.at((*itLS)->idnumberLS);
	if (_absDyn->isLdependent()){
          tmpamp *=_cachedDynIdLSMap.at((*itLS)->L).at(_absDyn->grandMaId(grandmaAmp));
        }
        amp+=tmpamp;
      }
      if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData);
      result+=amp;     
    }
  }

  if (!_absDyn->isLdependent()) result *=_cachedDynIdMap.at(_absDyn->grandMaId(grandmaAmp));
  result*=_isospinCG;
  return result;
}


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
    double theMag=std::abs(fitPar->Value(magName));
    double thePhi=fitPar->Value(phiName);

    _currentParamMagExpi[*itLS]=std::polar(theMag, thePhi);    
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
   _cachedDynIdLSMap[(*it)->L][_absDyn->grandMaId(grandmaAmp)]=_absDyn->eval(theData, grandmaAmp, (*it)->L);
 }  

 if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
 return;
}



