#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Examples/D0ToKsPipPim/D0ToKPiSPiDec.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
// #include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimEventList.hh"
#include "PwaDynamics/FVectorPiPiS.hh"

D0ToKPiSPiDec::D0ToKPiSPiDec(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates) :
  AbsXdecAmp(name, hypVec, 0)
  ,_iso1Key("KPiSWaveIso1")
  ,_iso3Key("KPiSWaveIso3")
  ,_iso1Hyp(false)
  ,_iso3Hyp(false)
  ,_theStatesPtr(theStates)
{
  initialize();
}

D0ToKPiSPiDec::~D0ToKPiSPiDec()
{
}

complex<double> D0ToKPiSPiDec::XdecAmp(Spin lamX, EvtDataNew* theData){

  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate){
    return _cachedAmpMap[theData->evtNo][lamX];
  }
  
  
  complex<double> result(0.,0.);
  complex<double> KpiSFVecMass(0.,0.);
  Vector4<double > p4KPi=theData->FourVecsDec[enumD0KsPiPi4V::V4_KsPim_HeliD0];
  
#ifdef _OPENMP
#pragma omp critical
  {
#endif
    _fVec->evalMatrix(p4KPi.M());    
    KpiSFVecMass=(*_fVec)[0];
#ifdef _OPENMP
  }
#endif
  
  result = KpiSFVecMass*complex<double> (1.,0.);


  if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical
    {
#endif
      _cachedAmpMap[evtNo][lamX]=result;
#ifdef _OPENMP
    }
#endif
  }

  return result;

}

void D0ToKPiSPiDec::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){
  std::map<std::string, double>::const_iterator it;
    for(it=_currentbFactorMap.begin();it!=_currentbFactorMap.end(); ++it){ 
      fitVal.otherParams[it->first]=it->second;
      fitErr.otherParams[it->first]=1.0;
    }

    for(it=_currentaProdMap.begin();it!=_currentaProdMap.end(); ++it){ 
      fitVal.otherParams[it->first]=it->second;
      fitErr.otherParams[it->first]=1.0;
    }

    for(it=_currentbProdMap.begin();it!=_currentbProdMap.end(); ++it){ 
      fitVal.otherParams[it->first]=it->second;
      fitErr.otherParams[it->first]=1.0;
    }

    for(it=_currentcProdMap.begin();it!=_currentcProdMap.end(); ++it){ 
      fitVal.otherParams[it->first]=it->second;
      fitErr.otherParams[it->first]=1.0;
    }

    for(it=_currentphaseProdMap.begin();it!=_currentphaseProdMap.end(); ++it){ 
      fitVal.otherParams[it->first]=it->second;
      fitErr.otherParams[it->first]=0.3;
    }

}

void D0ToKPiSPiDec::print(std::ostream& os) const{
  return; //dummy
}

void D0ToKPiSPiDec::initialize(){


    _currentaProdMap[_name+"a_Kpi"]=1.;

    _currentbProdMap[_name+"b_Kpi"]=1.;

    _currentcProdMap[_name+"c_Kpi"]=1.;

    _currentphaseProdMap[_name+"phi_KpiPosNeg"]=0.;

  if(_name==_iso1Key){
    _iso1Hyp=true;
    _currentbFactorMap[_name+"b_pole1Mag"]=1.;
    _currentbFactorMap[_name+"b_pole1PhiPosNeg"]=0.;
    _currentaProdMap[_name+"a_Ketap"]=0.5;
    _currentbProdMap[_name+"b_Ketap"]=0.5;
    _currentcProdMap[_name+"c_Ketap"]=0.5;
    _currentphaseProdMap[_name+"phi_KetapPosNeg"]=0.;
    _kMatr=  boost::shared_ptr<KMatrixKPiSFocus> (new KMatrixKPiSFocus(1));
  }
  else if(_name==_iso3Key){
    _iso3Hyp=true;
    _kMatr=  boost::shared_ptr<KMatrixKPiSFocus> (new KMatrixKPiSFocus(3));
  }
  else{
    Alert << "KPiS wave doesn't exist for key\t" << _name <<endmsg;
    exit(0);
  }

  _pVec = boost::shared_ptr<PVectorKPiSFocus> (new PVectorKPiSFocus(_kMatr));
  _fVec =  boost::shared_ptr<FVector> (new FVector(_kMatr, _pVec));

  Info << "hypothesis\t" << _name << "\t found" << endmsg;
}

void D0ToKPiSPiDec::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=false;

 std::map<std::string, double>::const_iterator it;
 if(_iso1Hyp){
    for(it=_currentbFactorMap.begin();it!=_currentbFactorMap.end(); ++it){
      double currentbFactor=theParamVal.otherParams[it->first];
      if ( fabs(currentbFactor-it->second) > 1.e-10){
	DebugMsg << "bFactor " << it->first << ":\t" << "current: " << it->second << "\tnew: " << currentbFactor << endmsg;
	_recalculate=true;
      } 
    }
 }

    for(it=_currentaProdMap.begin();it!=_currentaProdMap.end(); ++it){
      double currentaProd=theParamVal.otherParams[it->first];
      if ( fabs(currentaProd-it->second) > 1.e-10){
	DebugMsg << "aProd " << it->first << ":\t" << "current: " << it->second << "\tnew: " << currentaProd << endmsg;
	_recalculate=true;
      } 
    }

    for(it=_currentbProdMap.begin();it!=_currentbProdMap.end(); ++it){
      double currentbProd=theParamVal.otherParams[it->first];
      if ( fabs(currentbProd-it->second) > 1.e-10){
	DebugMsg << "bProd " << it->first << ":\t" << "current: " << it->second << "\tnew: " << currentbProd << endmsg;
	_recalculate=true;
      } 
    }

    for(it=_currentcProdMap.begin();it!=_currentcProdMap.end(); ++it){
      double currentcProd=theParamVal.otherParams[it->first];
      if ( fabs(currentcProd-it->second) > 1.e-10){
	DebugMsg << "cProd " << it->first << ":\t" << "current: " << it->second << "\tnew: " << currentcProd << endmsg;
	_recalculate=true;
      } 
    }

    for(it=_currentphaseProdMap.begin();it!=_currentphaseProdMap.end(); ++it){
      double currentphaseProd=theParamVal.otherParams[it->first];
      if ( fabs(currentphaseProd-it->second) > 1.e-10){
	DebugMsg << "phaseProd " << it->first << ":\t" << "current: " << it->second << "\tnew: " << currentphaseProd << endmsg;
	_recalculate=true;
      } 
    }
  
  if (_recalculate) Info << "Recalculate amplitude:\t" << _name << endmsg;
  
  
}


void D0ToKPiSPiDec::updateFitParams(fitParamsNew& theParamVal){
 
  std::map<std::string, double>::iterator it;

  if(_iso1Hyp){
    for(it=_currentbFactorMap.begin();it!=_currentbFactorMap.end(); ++it){
      double currentbFactor=theParamVal.otherParams[it->first];
      it->second = currentbFactor;
    }
  
  
    complex<double> b_pole1=_currentbFactorMap[_name+"b_pole1Mag"]*complex<double>(cos(_currentbFactorMap[_name+"b_pole1PhiPosNeg"]), sin(_currentbFactorMap[_name+"b_pole1PhiPosNeg"]));
    _pVec->updateBeta(0, b_pole1);
  }

  for(it=_currentaProdMap.begin();it!=_currentaProdMap.end(); ++it){
    double currentaFactor=theParamVal.otherParams[it->first];
    it->second = currentaFactor;
  }
  
  for(it=_currentbProdMap.begin();it!=_currentbProdMap.end(); ++it){
     double currentbFactor=theParamVal.otherParams[it->first];
     it->second = currentbFactor;
   }

  for(it=_currentcProdMap.begin();it!=_currentcProdMap.end(); ++it){
     double currentcFactor=theParamVal.otherParams[it->first];
     it->second = currentcFactor;
   }

  for(it=_currentphaseProdMap.begin();it!=_currentphaseProdMap.end(); ++it){
     double currentphaseFactor=theParamVal.otherParams[it->first];
     it->second = currentphaseFactor;
   }

  _pVec->updateAprod(0,_currentaProdMap[_name+"a_Kpi"]);
  _pVec->updateBprod(0,_currentbProdMap[_name+"b_Kpi"]);
  _pVec->updateCprod(0,_currentcProdMap[_name+"c_Kpi"]);
  _pVec->updatePhaseprod(0,_currentphaseProdMap[_name+"phi_KpiPosNeg"]);
 
  if(_iso1Hyp){
  _pVec->updateAprod(1,_currentaProdMap[_name+"a_Ketap"]);
  _pVec->updateBprod(1,_currentbProdMap[_name+"b_Ketap"]);
  _pVec->updateCprod(1,_currentcProdMap[_name+"c_Ketap"]);
  _pVec->updatePhaseprod(1,_currentphaseProdMap[_name+"phi_KetapPosNeg"]);
  }
}
