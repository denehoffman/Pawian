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

// KPiSWaveIso32Dynamics class definition file. -*- C++ -*-
// Copyright 2013 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/KPiSWaveIso32Dynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PVectorKPiSFocus.hh"
#include "FitParams/AbsPawianParameters.hh"

KPiSWaveIso32Dynamics::KPiSWaveIso32Dynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
{
  _kMatr =  std::shared_ptr<KMatrixKPiSFocus> (new KMatrixKPiSFocus(3));
  _isLdependent=false;
}

KPiSWaveIso32Dynamics::~KPiSWaveIso32Dynamics()
{
}

complex<double> KPiSWaveIso32Dynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();

  if ( _cacheAmps && !_recalcMap.at(currentKey)){
    result=_cachedStringMap.at(evtNo).at(currentKey);
  }

  else{
      theMutex.lock();
      result=_fVecMap.at(currentKey)->evalProjMatrix(theData->DoubleMassId.at(_dynId), 0);
      if ( _cacheAmps){
	_cachedStringMap[evtNo][currentKey]=result;
      }
      theMutex.unlock();
  }
  
  return result;
}

void KPiSWaveIso32Dynamics::fillDefaultParams(std::shared_ptr<AbsPawianParameters> fitPar){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentaProdMap.begin(); it1!=_currentaProdMap.end(); ++it1){
    std::string theName=it1->first;    
    std::map<std::string, double>::iterator it2;
 
    // a prod factors pure real   
    std::map<std::string, double>& aProds=it1->second;
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

void KPiSWaveIso32Dynamics::fillParamNameList(){
  _paramNameList.clear();
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentaProdMap.begin(); it1!=_currentaProdMap.end(); ++it1){
    std::string theName=it1->first;
    std::vector<std::string> currentNameList;    
    
    std::map<std::string, double>::iterator it2;
 
    // a prod factors pure real   
    std::map<std::string, double>& aProds=it1->second;
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

bool KPiSWaveIso32Dynamics::checkRecalculation(std::shared_ptr<AbsPawianParameters> fitParNew, std::shared_ptr<AbsPawianParameters> fitParOld){

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

void KPiSWaveIso32Dynamics::updateFitParams(std::shared_ptr<AbsPawianParameters> fitPar){

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentaProdMap.begin(); it1!=_currentaProdMap.end(); ++it1){
    std::string theName=it1->first;    
    std::map<std::string, double>::iterator it2;
 
    // a prod factors pure real   
    std::map<std::string, double>& aProds=it1->second;
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
    currentPVec->updateAprod(0, aProds.at("a"));
    currentPVec->updateBprod(0, bProds.at("b"));
    currentPVec->updateCprod(0, cProds.at("c"));
    currentPVec->updatePhaseprod(0, phaseProds.at("Phi"));
  }
}

void KPiSWaveIso32Dynamics::addGrandMa(std::shared_ptr<AbsDecay> theDec){
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

  _currentaProdMap[theName]["a"]=1.;

  _currentbProdMap[theName]["b"]=0.5;
  
  _currentcProdMap[theName]["c"]=0.1;
  
  _currentphaseProdMap[theName]["Phi"]=0.;

  currentPVector->updateAprod(0, _currentaProdMap[theName]["a"]);
  currentPVector->updateBprod(0, _currentbProdMap[theName]["b"]);
  currentPVector->updateCprod(0, _currentcProdMap[theName]["c"]);
  currentPVector->updatePhaseprod(0, _currentphaseProdMap[theName]["Phi"]);

  std::shared_ptr<FVector> currentFVector=std::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& KPiSWaveIso32Dynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  return grandmaAmp->absDec()->massParKey();
}

const unsigned short KPiSWaveIso32Dynamics::grandMaId(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaId;
  return grandmaAmp->absDec()->massParId(); 
}
