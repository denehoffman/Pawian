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

KPiSWaveIso32Dynamics::KPiSWaveIso32Dynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
{
  _kMatr =  boost::shared_ptr<KMatrixKPiSFocus> (new KMatrixKPiSFocus(3));
}

KPiSWaveIso32Dynamics::~KPiSWaveIso32Dynamics()
{
}

complex<double> KPiSWaveIso32Dynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";

  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->massParKey();

  if ( _cacheAmps && !_recalcMap[currentKey]){
    result=_cachedStringMap[evtNo][currentKey];
  }

  else{
      theMutex.lock();
      boost::shared_ptr<FVector> currentFVec=_fVecMap[currentKey];
      currentFVec->evalMatrix(theData->FourVecsString[_dynKey].M());
      result=(*currentFVec)[0];
      if ( _cacheAmps){
	_cachedStringMap[evtNo][currentKey]=result;
      }
      theMutex.unlock();
  }
  
  return result;
}

void  KPiSWaveIso32Dynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentaProdMap.begin(); it1!=_currentaProdMap.end(); ++it1){
    
    std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& aProds=it1->second;
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.;
    }

    std::map<std::string, double>& bProds=_currentbProdMap[it1->first];
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.; 
    }

    std::map<std::string, double>& cProds=_currentcProdMap[it1->first];
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.;
    }

    std::map<std::string, double>& phaseProds=_currentphaseProdMap[it1->first];
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.;
    }
  }
}

bool KPiSWaveIso32Dynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentaProdMap.begin(); it1!=_currentaProdMap.end(); ++it1){

    _recalcMap[it1->first]=false;    
    
     std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& aProds=it1->second;
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      if (fabs(it2->second - theParamVal.otherParams[it1->first+it2->first]) > 1.e-10){
	_recalculate=true;
	_recalcMap[it1->first]=true;
      }
    }

    std::map<std::string, double>& bProds=_currentbProdMap[it1->first];
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      if (fabs(it2->second - theParamVal.otherParams[it1->first+it2->first]) > 1.e-10){
	_recalculate=true;
	_recalcMap[it1->first]=true;
      }
    }

    std::map<std::string, double>& cProds=_currentcProdMap[it1->first];
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      if (fabs(it2->second - theParamVal.otherParams[it1->first+it2->first]) > 1.e-10){
	_recalculate=true;
	_recalcMap[it1->first]=true;
      }
    }

    std::map<std::string, double>& phaseProds=_currentphaseProdMap[it1->first];
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      if (fabs(it2->second - theParamVal.otherParams[it1->first+it2->first]) > 1.e-10){
	_recalculate=true;
	_recalcMap[it1->first]=true;
      }
    }

  }
  return _recalculate;
}

void KPiSWaveIso32Dynamics::updateFitParams(fitParams& theParamVal){

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentaProdMap.begin(); it1!=_currentaProdMap.end(); ++it1){

    std::map<std::string, double>::iterator it2;

    std::map<std::string, double>& aProds = it1->second;
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it1->first+it2->first];
    }

    std::map<std::string, double>& bProds=_currentbProdMap[it1->first];
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it1->first+it2->first];
    }

    std::map<std::string, double>& cProds=_currentcProdMap[it1->first];
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it1->first+it2->first];
    }

    std::map<std::string, double>& phaseProds=_currentphaseProdMap[it1->first];
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it1->first+it2->first];
    }

    boost::shared_ptr<PVectorKPiSFocus> currentPVec=_pVecMap[it1->first];
    currentPVec->updateAprod(0, aProds["a_PosNeg"]);
    currentPVec->updateBprod(0, bProds["b_PosNeg"]);
    currentPVec->updateCprod(0, cProds["c_PosNeg"]);
    currentPVec->updatePhaseprod(0, phaseProds["Phi"]);
  }
}

void KPiSWaveIso32Dynamics::addGrandMa(boost::shared_ptr<AbsDecay> theDec){
  if(0==theDec){
    Alert << "Can not add AbsXdecAmp; 0 pointer!!!" << endmsg;
    exit(1);
  }
  
  std::string theName=_massKey+theDec->massParKey();

  std::cout << "addGrandMa:\t" << theName << std::endl;

  std::map<std::string, boost::shared_ptr<FVector> >::iterator it = _fVecMap.find(theName);
  
  if (it != _fVecMap.end()) return;

  boost::shared_ptr<PVectorKPiSFocus> currentPVector=boost::shared_ptr<PVectorKPiSFocus>(new PVectorKPiSFocus(_kMatr));
  _pVecMap[theName]=currentPVector;

  _currentaProdMap[theName]["a_PosNeg"]=1.;

  _currentbProdMap[theName]["b_PosNeg"]=0.5;
  
  _currentcProdMap[theName]["c_PosNeg"]=0.1;
  
  _currentphaseProdMap[theName]["Phi"]=0.;

  currentPVector->updateAprod(0, _currentaProdMap[theName]["a_PosNeg"]);
  currentPVector->updateBprod(0, _currentbProdMap[theName]["b_PosNeg"]);
  currentPVector->updateCprod(0, _currentcProdMap[theName]["c_PosNeg"]);
  currentPVector->updatePhaseprod(0, _currentphaseProdMap[theName]["Phi"]);

  boost::shared_ptr<FVector> currentFVector=boost::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& KPiSWaveIso32Dynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  return grandmaAmp->absDec()->massParKey();
}
