// KPiSWaveDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/KPiSWaveDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Particle/Particle.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/PVectorKPiSFocus.hh"

KPiSWaveDynamics::KPiSWaveDynamics(std::string& name, std::vector<Particle*>& fsParticles, Particle* mother) :
  AbsDynamics(name, fsParticles, mother)
{
  _kMatr =  boost::shared_ptr<KMatrixKPiSFocus> (new KMatrixKPiSFocus(1));
}

KPiSWaveDynamics::~KPiSWaveDynamics()
{
}

complex<double> KPiSWaveDynamics::eval(EvtData* theData, AbsXdecAmp* grandmaAmp, Spin OrbMom){

complex<double> result(0.,0.);

  int evtNo=theData->evtNo;
  std::string currentKey="default";
  if(0!=grandmaAmp) currentKey=_massKey+grandmaAmp->absDec()->dynKey(); 

  if ( _cacheAmps && !_recalcMap[currentKey]){
    result=_cachedStringMap[evtNo][currentKey];
  }

  else{
#ifdef _OPENMP
#pragma omp critical (KPiSWaveDynCache)
    {
#endif
      boost::shared_ptr<FVector> currentFVec=_fVecMap[currentKey];
      currentFVec->evalMatrix(theData->FourVecsString[_dynKey].M());
      result=(*currentFVec)[0];
      if ( _cacheAmps){
	_cachedStringMap[evtNo][currentKey]=result;
      }
#ifdef _OPENMP
    }
#endif
  }
  
  return result;
}

void  KPiSWaveDynamics::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    
    std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      fitVal.otherParams[it1->first+it2->first]=it2->second;
      fitErr.otherParams[it1->first+it2->first]=1.;
    }

    std::map<std::string, double>& aProds=_currentaProdMap[it1->first];
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

bool KPiSWaveDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){

    _recalcMap[it1->first]=false;    
    
     std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double>& bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      if (fabs(it2->second - theParamVal.otherParams[it1->first+it2->first]) > 1.e-10){
	_recalculate=true;
	_recalcMap[it1->first]=true;
      }
    }

    std::map<std::string, double>& aProds=_currentaProdMap[it1->first];
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

void KPiSWaveDynamics::updateFitParams(fitParams& theParamVal){

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){

    std::map<std::string, double>::iterator it2;

    it1->second["b_pole1Mag"]=theParamVal.otherParams[it1->first+"b_pole1Mag"];
    it1->second["b_pole1Phi"]=theParamVal.otherParams[it1->first+"b_pole1Phi"];

    complex<double> b_pole1=it1->second["b_pole1Mag"]*complex<double>(cos(it1->second["b_pole1Phi"]), sin(it1->second["b_pole1Phi"]));

    std::map<std::string, double>& aProds=_currentaProdMap[it1->first];
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
    currentPVec->updateBeta(0, b_pole1);
    currentPVec->updateAprod(0, aProds["a_KpiPosNeg"]);
    currentPVec->updateBprod(0, bProds["b_KpiPosNeg"]);
    currentPVec->updateCprod(0, cProds["c_KpiPosNeg"]);
    currentPVec->updatePhaseprod(0, phaseProds["KpiPhi"]);
    currentPVec->updateAprod(1, aProds["a_KetapPosNeg"]);
    currentPVec->updateBprod(1, bProds["b_KetapPosNeg"]);
    currentPVec->updateCprod(1, cProds["c_KetapPosNeg"]);
    currentPVec->updatePhaseprod(1, phaseProds["KetapPhi"]);

  }
}

void KPiSWaveDynamics::addGrandMa(boost::shared_ptr<AbsDecay> theDec){
  if(0==theDec){
    Alert << "Can not add AbsXdecAmp; 0 pointer!!!" << endmsg;
    exit(1);
  }
  
  //  std::string theName=_massKey+theDec->motherJPC()->name();
  std::string theName=_massKey+theDec->dynKey();
  std::cout << "addGrandMa:\t" << theName << std::endl;

  std::map<std::string, boost::shared_ptr<FVector> >::iterator it = _fVecMap.find(theName);
  
  if (it != _fVecMap.end()) return;

  boost::shared_ptr<PVectorKPiSFocus> currentPVector=boost::shared_ptr<PVectorKPiSFocus>(new PVectorKPiSFocus(_kMatr));
  _pVecMap[theName]=currentPVector;

  _currentaProdMap[theName]["a_KpiPosNeg"]=1.;
  _currentaProdMap[theName]["a_KetapPosNeg"]=1.;  

  _currentbProdMap[theName]["b_KpiPosNeg"]=1.;
  _currentbProdMap[theName]["b_KetapPosNeg"]=0.5;
  
  _currentcProdMap[theName]["c_KpiPosNeg"]=1.;
  _currentcProdMap[theName]["c_KetapPosNeg"]=0.5;  
  
  _currentphaseProdMap[theName]["KpiPhi"]=0.;
  _currentphaseProdMap[theName]["KetapPhi"]=0.;

  _currentbFactorMap[theName]["b_pole1Mag"]=1.;
  _currentbFactorMap[theName]["b_pole1Phi"]=0.;

  complex<double> b_pole1=_currentbFactorMap[theName]["b_pole1Mag"]*complex<double>(cos(_currentbFactorMap[theName]["b_pole1Phi"]), sin(_currentbFactorMap[theName]["b_pole1Phi"]));
  currentPVector->updateBeta(0, b_pole1);
  currentPVector->updateAprod(0, _currentaProdMap[theName]["a_KpiPosNeg"]);
  currentPVector->updateBprod(0, _currentbProdMap[theName]["b_KpiPosNeg"]);
  currentPVector->updateCprod(0, _currentcProdMap[theName]["c_KpiPosNeg"]);
  currentPVector->updatePhaseprod(0, _currentphaseProdMap[theName]["KpiPhi"]);
  currentPVector->updateAprod(1, _currentaProdMap[theName]["a_KetapPosNeg"]);
  currentPVector->updateBprod(1, _currentbProdMap[theName]["b_KetapPosNeg"]);
  currentPVector->updateCprod(1, _currentcProdMap[theName]["c_KetapPosNeg"]);
  currentPVector->updatePhaseprod(1, _currentphaseProdMap[theName]["KetapPhi"]);

  boost::shared_ptr<FVector> currentFVector=boost::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;
  _recalcMap[theName]=true;
}

const std::string& KPiSWaveDynamics::grandMaKey(AbsXdecAmp* grandmaAmp){
  if(0==grandmaAmp) return _grandmaKey;
  return grandmaAmp->absDec()->dynKey();
}
