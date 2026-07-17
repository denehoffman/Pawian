//************************************************************************//
//									  //
//  Copyright 2014 Bertram Kopf (bertram@ep1.rub.de)			  //
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

// LSDecAmps class definition file. -*- C++ -*-
// Copyright 2014 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/LSDecAmps.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "PwaUtils/DataUtils.hh"
#include "PwaUtils/GlobalEnv.hh"
#include "PwaUtils/IsobarLSDecay.hh"
#include "PwaDynamics/BarrierFactor.hh"
#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"
#include "ErrLogger/ErrLogger.hh"
#include "FitParams/AbsPawianParameters.hh"

LSDecAmps::LSDecAmps(std::shared_ptr<IsobarLSDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_LSs(theDec->LSAmps())
  ,_factorMag(1.)
   ,_Smax(0)
{
  initialize();
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    if( (*it)->S > _Smax ) _Smax=(*it)->S;
  }

  if(_LSs.size()>0) _factorMag=1./sqrt(_LSs.size());
  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  InfoMsg << "_parityFactor=\t" << _parityFactor << endmsg;
  fillCgPreFactor();
}

LSDecAmps::LSDecAmps(std::shared_ptr<AbsDecay> theDec, ChannelID channelID) :
  AbsXdecAmp(theDec, channelID)
  ,_Smax(0)
{
  initialize();
  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    if( (*it)->S > _Smax ) _Smax=(*it)->S;
  }

  Particle* daughter1=_decay->daughter1Part();
  Particle* daughter2=_decay->daughter2Part();
  _parityFactor=daughter1->theParity()*daughter2->theParity()*pow(-1,_JPCPtr->J-daughter1->J()-daughter2->J());
  InfoMsg << "_parityFactor=\t" << _parityFactor << endmsg;
  fillCgPreFactor();
}

LSDecAmps::~LSDecAmps()
{
}

complex<double> LSDecAmps::XdecPartAmp(const Spin& lamX, Spin& lamDec, short fixDaughterNr, 
				       EvtData* theData, Spin& lamFs, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
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

   result=lsLoop( grandmaAmp, lamX, theData, lam1Min, lam1Max, lam2Min, lam2Max, false);

  return result;
}

complex<double> LSDecAmps::productionHelicityComponent(
    const Spin& totalProjection, const Spin& producedProjection,
    const Spin& recoilHelicity, EvtData* theData, AbsXdecAmp* grandmaAmp) {
  Spin producedMin = producedProjection;
  Spin producedMax = producedProjection;
  Spin recoilMin = recoilHelicity;
  Spin recoilMax = recoilHelicity;
  return lsLoop(grandmaAmp, totalProjection, theData, producedMin, producedMax,
                recoilMin, recoilMax, false);
}

complex<double> LSDecAmps::XdecAmp(const Spin& lamX, EvtData* theData, AbsXdecAmp* grandmaAmp){

  complex<double> result(0.,0.);
  if( fabs(lamX) > _JPCPtr->J) return result;

  short currentSpinIndex=FunctionUtils::spin1IdIndex(_projId, lamX);

  if (!_recalculate){
    return _cachedAmpIdMap.at(theData->evtNo).at(_absDyn->grandMaId(grandmaAmp)).at(currentSpinIndex);
  }

  result=lsLoop(grandmaAmp, lamX, theData, _lam1MinProj, _lam1MaxProj, _lam2MinProj, _lam2MaxProj, true);
 
  if ( _cacheAmps){
      _cachedAmpIdMap[theData->evtNo][_absDyn->grandMaId(grandmaAmp)][currentSpinIndex]=result;
   }

  if(result.real()!=result.real()){
    InfoMsg << "dyn name: " << _absDyn->name() 
	    << "\nname(): " << name()
	    << endmsg;
    Alert << "result:\t" << result << endmsg;
    exit(0);
  }
  return result;
}

complex<double> LSDecAmps::lsLoop(AbsXdecAmp* grandmaAmp, const Spin& lamX, EvtData* theData, 
				  Spin& lam1Min, Spin& lam1Max, Spin& lam2Min, Spin& lam2Max, bool withDecs){
 
  complex<double> result(0.,0.);

  std::map<Id3StringType, complex<double> >& currentWignerDMap=theData->WignerDIdId3.at(_decay->wigDWigDRefId());
  std::map< Spin, std::map<Spin, std::map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > > >::iterator itLamLamLS;
  std::map<Spin, std::map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > >::iterator itLamLS;
  std::map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >::iterator itLS;

  for(itLamLamLS=_cgPreFactor_LamLamLSMap.begin(); itLamLamLS!=_cgPreFactor_LamLamLSMap.end(); ++itLamLamLS){
    //Spin lambda1(itLamLamLS->first);
    const Spin& lambda1=itLamLamLS->first; 
    if(lambda1<lam1Min || lambda1>lam1Max) continue;
    std::map<Spin, std::map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess > >& cgPre_lam2LSMap = itLamLamLS->second;
    for(itLamLS=cgPre_lam2LSMap.begin(); itLamLS!=cgPre_lam2LSMap.end(); ++itLamLS){
      //Spin lambda2(itLamLS->first);
      const Spin& lambda2=itLamLS->first;
      if(lambda2<lam2Min || lambda2>lam2Max) continue;
      Spin lambda(lambda1-lambda2);                                                                  
      // if( std::abs(lambda)>_JPCPtr->J || std::abs(lambda)>_Smax) continue; already done in LSDecAmps::fillCgPreFactor()
      complex<double> amp(0.,0.);
      std::map<std::shared_ptr<const LScomb>, double, pawian::Collection::SharedPtrLess >& cgPre_LSMap = itLamLS->second;
      for(itLS=cgPre_LSMap.begin(); itLS!=cgPre_LSMap.end(); ++itLS){
	// if( std::abs(lambda)> itLS->first->S) continue; already done in LSDecAmps::fillCgPreFactor()
        if (_absDyn->isLdependent()){
          amp += _currentParamPreFacMagExpi.at(itLS->first)*itLS->second
            * _cachedDynIdLSMap.at(itLS->first->L).at(_absDyn->grandMaId(grandmaAmp));
        }
        else amp+=_currentParamPreFacMagExpi.at(itLS->first)*itLS->second;
      }

      Id3StringType IdJLamXLam12=FunctionUtils::spin3Index(_J, lamX, lambda);
      amp *= conj(currentWignerDMap.at(IdJLamXLam12));
      if(withDecs) amp *=daughterAmp(lambda1, lambda2, theData);
      result+=amp;
    }
  }
  if (!_absDyn->isLdependent()){
    result *=_cachedDynIdMap.at(_absDyn->grandMaId(grandmaAmp));
  }
  return result;
}

void LSDecAmps::calcDynamics(EvtData* theData, AbsXdecAmp* grandmaAmp){
  if(!_recalculate) return; 

  if(!_absDyn->isLdependent()){
    AbsXdecAmp::calcDynamics(theData, grandmaAmp);
    //    return;
  }
  else{
    std::vector< std::shared_ptr<const LScomb> >::iterator it;
    for (it=_LSs.begin(); it!=_LSs.end(); ++it){
      _cachedDynIdLSMap[(*it)->L][_absDyn->grandMaId(grandmaAmp)]
	= _absDyn->eval(theData, grandmaAmp, (*it)->L);
    }
  }

 if(!_daughter1IsStable) _decAmpDaughter1->calcDynamics(theData, this);
 if(!_daughter2IsStable) _decAmpDaughter2->calcDynamics(theData, this);
 return;
}

void  LSDecAmps::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double valMag=_factorMag;
    double errMag=_factorMag/2.;

    fitPar->Add(magName, valMag, errMag);

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

void LSDecAmps::fillParamNameList(){
  _paramNameList.clear();
  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    _paramNameList.push_back(magName);
    //phi
    std::string phiName=(*itLS)->name()+_key+"Phi";
    _paramNameList.push_back(phiName);
  }

}

void LSDecAmps::print(std::ostream& os) const{
  return; //dummy
}


void LSDecAmps::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::vector< std::shared_ptr<const LScomb> >::const_iterator itLS;
  for(itLS=_LSs.begin(); itLS!=_LSs.end(); ++itLS){
    //fill magnitude
    std::string magName=(*itLS)->name()+_key+"Mag";
    double theMag= std::abs(fitPar->Value(magName));
    
    std::string phiName=(*itLS)->name()+_key+"Phi";
    double thePhi=fitPar->Value(phiName);

    _currentParamPreFacMagExpi[*itLS]=_preFactor*_isospinCG*std::polar(theMag, thePhi);
  }

   _absDyn->updateFitParams(fitPar);

  if(!_daughter1IsStable) _decAmpDaughter1->updateFitParams(fitPar);
  if(!_daughter2IsStable) _decAmpDaughter2->updateFitParams(fitPar);
}

void  LSDecAmps::fillCgPreFactor(){

  std::vector< std::shared_ptr<const LScomb> >::iterator it;
  for (it=_LSs.begin(); it!=_LSs.end(); ++it){
    for(Spin lambda1=-_Jdaughter1; lambda1<=_Jdaughter1; ++lambda1){
      for(Spin lambda2=-_Jdaughter2; lambda2<=_Jdaughter2; ++lambda2){
	Spin lambda = lambda1-lambda2;
	if( fabs(lambda)>_JPCPtr->J || fabs(lambda)>(*it)->S) continue;
	if( std::abs(lambda)> (*it)->S) continue;
	double currentResult=sqrt(2.*(*it)->L+1)
          *Clebsch((*it)->L, 0, (*it)->S, lambda, _JPCPtr->J, lambda)
          *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  );

	if(norm(currentResult)>1.e-10) _cgPreFactor_LamLamLSMap[lambda1][lambda2][*it]=currentResult;
	
	//	_cgPreFactor_LamLamLSMap[lambda1][lambda2][*it]=sqrt(2.*(*it)->L+1)
	//	  *Clebsch((*it)->L, 0, (*it)->S, lambda, _JPCPtr->J, lambda)
	//	  *Clebsch(_Jdaughter1, lambda1, _Jdaughter2, -lambda2, (*it)->S, lambda  );
      }
    }
  }
}
