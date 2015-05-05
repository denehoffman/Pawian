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

// KPiSWaveIso12Dynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/KPiSWaveIso12Dynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PVectorKPiSFocus.hh"
#include "FitParams/AbsPawianParameters.hh"

KPiSWaveIso12Dynamics::KPiSWaveIso12Dynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
{
  _kMatr =  std::shared_ptr<KMatrixKPiSFocus> (new KMatrixKPiSFocus(1));
  _isLdependent=false;
}

KPiSWaveIso12Dynamics::~KPiSWaveIso12Dynamics()
{
}

complex<double> KPiSWaveIso12Dynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();

  if ( _cacheAmps && !_recalcMap.at(currentKey)){
    result=_cachedStringMap.at(evtNo).at(currentKey);
  }

  else{
      theMutex.lock();
      result=_fVecMap.at(currentKey)->evalProjMatrix(theData->DoubleString.at(_dynKey), 0);
      if ( _cacheAmps){
	_cachedStringMap[evtNo][currentKey]=result;
      }
      theMutex.unlock();
  }
  
  return result;
}


void KPiSWaveIso12Dynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMagMap.begin(); it1!=_currentbFactorMagMap.end(); ++it1){
    std::string theName=it1->first;
    std::map<std::string, double>& bMagFactors = it1->second;
    
    std::map<std::string, double>::iterator it2;
    for(it2=bMagFactors.begin(); it2!=bMagFactors.end(); ++it2){
      std::string magName=theName+it2->first;
      fitPar->Add(magName, it2->second , 1.);
      fitPar->SetLimits(magName, 0., it2->second+30.);
    }

    std::map<std::string, double>& bPhiFactors = _currentbFactorPhiMap.at(it1->first);
    for(it2=bPhiFactors.begin(); it2!=bPhiFactors.end(); ++it2){
      std::string phiName=theName+it2->first;
      fitPar->Add(phiName, 0. , 0.2);
    }

    // a prod factors pure real
    std::map<std::string, double>& aProds=_currentaProdMap.at(it1->first);
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      std::string aProdName=theName+it2->first;
      fitPar->Add(aProdName, it2->second , 0.5);    
    }

    // b prod factors pure real
    std::map<std::string, double>& bProds=_currentbProdMap.at(it1->first);
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      std::string bProdName=theName+it2->first;
      fitPar->Add(bProdName, it2->second , 0.5);    
    }

    // c prod factors pure real
    std::map<std::string, double>& cProds=_currentcProdMap.at(it1->first);
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      std::string cProdName=theName+it2->first;
      fitPar->Add(cProdName, it2->second , 0.5);    
    }

    // phase prod factors
    std::map<std::string, double>& phaseProds=_currentphaseProdMap.at(it1->first);
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      std::string phaseName=theName+it2->first;
      fitPar->Add(phaseName, it2->second , 0.2);    
    }
  }
}
void KPiSWaveIso12Dynamics::fillParamNameList(){

  _paramNameList.clear();

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMagMap.begin(); it1!=_currentbFactorMagMap.end(); ++it1){
    std::string theName=it1->first;

    std::vector<std::string> currentNameList;
    
    std::map<std::string, double>& bMagFactors = it1->second;
    
    std::map<std::string, double>::iterator it2;
    for(it2=bMagFactors.begin(); it2!=bMagFactors.end(); ++it2){
      std::string magName=theName+it2->first;
      _paramNameList.push_back(magName);
      currentNameList.push_back(magName);
    }

    std::map<std::string, double>& bPhiFactors = _currentbFactorPhiMap.at(it1->first);
    for(it2=bPhiFactors.begin(); it2!=bPhiFactors.end(); ++it2){
      std::string phiName=theName+it2->first;
      _paramNameList.push_back(phiName);
      currentNameList.push_back(phiName);
    }

    // a prod factors pure real
    std::map<std::string, double>& aProds=_currentaProdMap.at(it1->first);
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      std::string aProdName=theName+it2->first;
      _paramNameList.push_back(aProdName);
      currentNameList.push_back(aProdName);
    }

    // b prod factors pure real
    std::map<std::string, double>& bProds=_currentbProdMap.at(it1->first);
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      std::string bProdName=theName+it2->first;
      _paramNameList.push_back(bProdName);
      currentNameList.push_back(bProdName);
    } 
   // c prod factors pure real
    std::map<std::string, double>& cProds=_currentcProdMap.at(it1->first);
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      std::string cProdName=theName+it2->first;
      _paramNameList.push_back(cProdName);
      currentNameList.push_back(cProdName);
    }

    // phase prod factors
    std::map<std::string, double>& phaseProds=_currentphaseProdMap.at(it1->first);
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      std::string phaseName=theName+it2->first;
      _paramNameList.push_back(phaseName);
      currentNameList.push_back(phaseName);
    }
    _paramNameListMap[theName]=currentNameList;
  }

}

bool KPiSWaveIso12Dynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){

  std::map<std::string, std::vector<std::string> >::iterator itMap;
  for(itMap= _paramNameListMap.begin(); itMap !=_paramNameListMap.end(); ++itMap){
    _recalcMap.at(itMap->first)=false;
    std::vector<std::string>::iterator itStr;
    for (itStr=itMap->second.begin(); itStr!=itMap->second.end(); ++itStr){
      std::string currentParamName=*itStr;
      if(!CheckDoubleEquality(fitParNew->Value(currentParamName), fitParOld->Value(currentParamName))){
        _recalcMap.at(itMap->first)=true;
        continue;
      }
    }
  }

  return AbsParamHandler::checkRecalculation(fitParNew, fitParOld);
}


void KPiSWaveIso12Dynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMagMap.begin(); it1!=_currentbFactorMagMap.end(); ++it1){
    std::string theName=it1->first;

    std::string bmagName=theName+"b_pole1Mag";
    std::string bphiName=theName+"b_pole1Phi";

    it1->second["b_pole1Mag"]=fitPar->Value(bmagName);

    double currentbFactorPhi=fitPar->Value(bphiName);
    _currentbFactorPhiMap.at(it1->first).at("b_pole1Phi")=currentbFactorPhi;

    complex<double> b_pole1=it1->second["b_pole1Mag"]*complex<double>(cos(currentbFactorPhi), sin(currentbFactorPhi));    

    std::map<std::string, double>::iterator it2;
    // a prod factors pure real
    std::map<std::string, double>& aProds=_currentaProdMap.at(it1->first);
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      std::string aProdName=theName+it2->first;
      it2->second = fitPar->Value(aProdName);
    }

    // b prod factors pure real
    std::map<std::string, double>& bProds=_currentbProdMap.at(it1->first);
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      std::string bProdName=theName+it2->first;
      it2->second = fitPar->Value(bProdName);
    } 
   // c prod factors pure real
    std::map<std::string, double>& cProds=_currentcProdMap.at(it1->first);
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      std::string cProdName=theName+it2->first;
      it2->second = fitPar->Value(cProdName);
    }
    // phase prod factors
    std::map<std::string, double>& phaseProds=_currentphaseProdMap.at(it1->first);
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      std::string phaseName=theName+it2->first;
      it2->second = fitPar->Value(phaseName);    
    }

    std::shared_ptr<PVectorKPiSFocus> currentPVec=_pVecMap.at(it1->first);
    currentPVec->updateBeta(0, b_pole1);
    currentPVec->updateAprod(0, aProds["a_Kpi"]);
    currentPVec->updateBprod(0, bProds["b_Kpi"]);
    currentPVec->updateCprod(0, cProds["c_Kpi"]);
    currentPVec->updatePhaseprod(0, phaseProds["KpiPhi"]);
    currentPVec->updateAprod(1, aProds["a_Ketap"]);
    currentPVec->updateBprod(1, bProds["b_Ketap"]);
    currentPVec->updateCprod(1, cProds["c_Ketap"]);
    currentPVec->updatePhaseprod(1, phaseProds["KetapPhi"]);
  }
}



void KPiSWaveIso12Dynamics::addGrandMa(std::shared_ptr<AbsDecay> theDec){
  if(0==theDec){
    Alert << "Can not add AbsXdecAmp; 0 pointer!!!" << endmsg;
    exit(1);
  }
  
  std::string theName=_massKey+theDec->massParKey();

  std::cout << "addGrandMa:\t" << theName << std::endl;

  std::map<std::string, std::shared_ptr<FVector> >::iterator it = _fVecMap.find(theName);
  
  if (it != _fVecMap.end()) return;

  std::shared_ptr<PVectorKPiSFocus> currentPVector=std::shared_ptr<PVectorKPiSFocus>(new PVectorKPiSFocus(_kMatr));
  _pVecMap[theName]=currentPVector;

  _currentaProdMap[theName]["a_Kpi"]=1.;
  _currentaProdMap[theName]["a_Ketap"]=1.;  

  _currentbProdMap[theName]["b_Kpi"]=0.5;
  _currentbProdMap[theName]["b_Ketap"]=0.5;
  
  _currentcProdMap[theName]["c_Kpi"]=0.1;
  _currentcProdMap[theName]["c_Ketap"]=0.1;  
  
  _currentphaseProdMap[theName]["KpiPhi"]=0.;
  _currentphaseProdMap[theName]["KetapPhi"]=0.;

  _currentbFactorMagMap[theName]["b_pole1Mag"]=1.;
  _currentbFactorPhiMap[theName]["b_pole1Phi"]=0.;

  complex<double> b_pole1=_currentbFactorMagMap[theName]["b_pole1Mag"]*complex<double>(cos(_currentbFactorPhiMap[theName]["b_pole1Phi"]), sin(_currentbFactorPhiMap[theName]["b_pole1Phi"]));
  currentPVector->updateBeta(0, b_pole1);
  currentPVector->updateAprod(0, _currentaProdMap[theName]["a_Kpi"]);
  currentPVector->updateBprod(0, _currentbProdMap[theName]["b_Kpi"]);
  currentPVector->updateCprod(0, _currentcProdMap[theName]["c_Kpi"]);
  currentPVector->updatePhaseprod(0, _currentphaseProdMap[theName]["KpiPhi"]);
  currentPVector->updateAprod(1, _currentaProdMap[theName]["a_Ketap"]);
  currentPVector->updateBprod(1, _currentbProdMap[theName]["b_Ketap"]);
  currentPVector->updateCprod(1, _currentcProdMap[theName]["c_Ketap"]);
  currentPVector->updatePhaseprod(1, _currentphaseProdMap[theName]["KetapPhi"]);

  std::shared_ptr<FVector> currentFVector=std::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& KPiSWaveIso12Dynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  //  return grandmaAmp->absDec()->dynKey();
  return grandmaAmp->absDec()->massParKey();
}
