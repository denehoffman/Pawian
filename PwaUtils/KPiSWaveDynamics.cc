// KPiSWaveDynamics class definition file. -*- C++ -*-
// Copyright 20123Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/KPiSWaveDynamics.hh"
#include "PwaUtils/XdecAmpRegistry.hh"
#include "PwaUtils/AbsDecay.hh"
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

complex<double> KPiSWaveDynamics::eval(EvtData* theData, Spin OrbMom){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedMap[evtNo];
  }
  complex<double> result(0.,0.);
  
  if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical (KPiSWaveDynCache)
    {
#endif
      _cachedMap[evtNo]=result;
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
    
    std::map<std::string, double> bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      fitVal.otherParams[it2->first]=it2->second;
      fitErr.otherParams[it2->first]=1.;
    }

    std::map<std::string, double> aProds=_currentaProdMap[it1->first];
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      fitVal.otherParams[it2->first]=it2->second;
      fitErr.otherParams[it2->first]=1.;
    }

    std::map<std::string, double> bProds=_currentbProdMap[it1->first];
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      fitVal.otherParams[it2->first]=it2->second;
      fitErr.otherParams[it2->first]=1.; 
    }

    std::map<std::string, double> cProds=_currentcProdMap[it1->first];
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      fitVal.otherParams[it2->first]=it2->second;
      fitErr.otherParams[it2->first]=1.;
    }

    std::map<std::string, double> phaseProds=_currentphaseProdMap[it1->first];
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      fitVal.otherParams[it2->first]=it2->second;
      fitErr.otherParams[it2->first]=1.;
    }
  }
}

bool KPiSWaveDynamics::checkRecalculation(fitParams& theParamVal){
  _recalculate=false;

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){
    
    std::map<std::string, double>::iterator it2;
    
    std::map<std::string, double> bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      if (fabs(it2->second = theParamVal.otherParams[it2->first]) > 1.e-10){
	_recalculate=true;
	return _recalculate;
      }
    }

    std::map<std::string, double> aProds=_currentaProdMap[it1->first];
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      if (fabs(it2->second = theParamVal.otherParams[it2->first]) > 1.e-10){
	_recalculate=true;
	return _recalculate;
      }
    }

    std::map<std::string, double> bProds=_currentbProdMap[it1->first];
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      if (fabs(it2->second = theParamVal.otherParams[it2->first]) > 1.e-10){
	_recalculate=true;
	return _recalculate;
      }
    }

    std::map<std::string, double> cProds=_currentcProdMap[it1->first];
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      if (fabs(it2->second = theParamVal.otherParams[it2->first]) > 1.e-10){
	_recalculate=true;
	return _recalculate;
      }
    }

    std::map<std::string, double> phaseProds=_currentphaseProdMap[it1->first];
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      if (fabs(it2->second = theParamVal.otherParams[it2->first]) > 1.e-10){
	_recalculate=true;
	return _recalculate;
      }
    }

  }
  return _recalculate;
}

void KPiSWaveDynamics::updateFitParams(fitParams& theParamVal){

  std::map<std::string, std::map<std::string, double> >::iterator it1;
  for(it1=_currentbFactorMap.begin(); it1!=_currentbFactorMap.end(); ++it1){

    std::map<std::string, double>::iterator it2;

    std::map<std::string, double> bFactors=it1->second;
    for(it2=bFactors.begin(); it2!=bFactors.end(); ++it2){
      it2->second = theParamVal.otherParams[it2->first];
    }

    std::map<std::string, double> aProds=_currentaProdMap[it1->first];
    for(it2=aProds.begin(); it2!=aProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it2->first];
    }

    std::map<std::string, double> bProds=_currentbProdMap[it1->first];
    for(it2=bProds.begin(); it2!=bProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it2->first];
    }

    std::map<std::string, double> cProds=_currentcProdMap[it1->first];
    for(it2=cProds.begin(); it2!=cProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it2->first];
    }

    std::map<std::string, double> phaseProds=_currentphaseProdMap[it1->first];
    for(it2=phaseProds.begin(); it2!=phaseProds.end(); ++it2){
      it2->second = theParamVal.otherParams[it2->first];
    }
  }
}

void KPiSWaveDynamics::addGrandMa(boost::shared_ptr<AbsDecay> theDec){
  if(0==theDec){
    Alert << "Can not add AbsXdecAmp; 0 pointer!!!" << endmsg;
    exit(1);
  }
  
//  std::string theName=theDec->name();
//  std::string theName=_massKey+theDec->motherJPC()->name();
  std::string theName=_name+theDec->motherJPC()->name();
  std::cout << "addGrandMa:\t" << theName << std::endl;

//  std::map<std::string, boost::shared_ptr<AbsXdecAmp> >::iterator it = _grandMaAmpMap.find(theName);

//  if (it !=_grandMaAmpMap.end()) return; //already there

//  boost::shared_ptr<AbsXdecAmp> currentAmp=XdecAmpRegistry::instance()->getXdecAmp(theDec);
//  _grandMaAmpMap[theName]=currentAmp;

  std::map<std::string, boost::shared_ptr<FVector> >::iterator it = _fVecMap.find(theName);
  
  if (it != _fVecMap.end()) return;

  std::string currentKey=_massKey+theDec->motherJPC()->name();
  
  boost::shared_ptr<PVectorKPiSFocus> currentPVector=boost::shared_ptr<PVectorKPiSFocus>(new PVectorKPiSFocus(_kMatr));
  _pVecMap[theName]=currentPVector;
  boost::shared_ptr<FVector> currentFVector=boost::shared_ptr<FVector>(new FVector(_kMatr, currentPVector));
  _fVecMap[theName]=currentFVector;

  _currentaProdMap[theName][currentKey+"a_KpiPosNeg"]=1.;
  _currentaProdMap[theName][currentKey+"a_KetapPosNeg"]=1.;  

  _currentbProdMap[theName][currentKey+"b_KpiPosNeg"]=1.;
  _currentbProdMap[theName][currentKey+"b_KetapPosNeg"]=0.5;
  
  _currentcProdMap[theName][currentKey+"c_KpiPosNeg"]=1.;
  _currentcProdMap[theName][currentKey+"c_KetapPosNeg"]=0.5;  
  
  _currentphaseProdMap[theName][currentKey+"_KpiPhi"]=0.;
  _currentphaseProdMap[theName][currentKey+"_KetapPhi"]=0.;

  _currentbFactorMap[theName][currentKey+"b_pole1Mag"]=1.;
  _currentbFactorMap[theName][currentKey+"b_pole1Phi"]=0.;
}
