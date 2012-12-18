#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Examples/D0ToKsPipPim/D0ToPiPiSKDec.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
// #include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimEventList.hh"
#include "PwaDynamics/FVectorPiPiS.hh"

D0ToPiPiSKDec::D0ToPiPiSKDec(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates) :
  AbsXdecAmp(name, hypVec, 0)
  ,_ASKey("PiPiSWaveAS")
  ,_piPiSASHyp(false)
  ,_theStatesPtr(theStates)
  ,_pipiSFVec(new FVectorPiPiS())
  ,_currentS0Val(0.)
{
  initialize();
}

D0ToPiPiSKDec::~D0ToPiPiSKDec()
{
}

complex<double> D0ToPiPiSKDec::XdecAmp(Spin lamX, EvtDataNew* theData){

  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate) return _cachedAmpMap[theData->evtNo][lamX];
  
  complex<double> result(0.,0.);
  complex<double> pipiSFVecMass(0.,0.);
  Vector4<double > p4PiPi=theData->FourVecsDec[enumD0KsPiPi4V::V4_PipPim_HeliD0];
  
#ifdef _OPENMP
#pragma omp critical
  {
#endif
    _pipiSFVec->evalMatrix(p4PiPi.M());    
    pipiSFVecMass=(*_pipiSFVec)[0];
#ifdef _OPENMP
  }
#endif
  
  result = pipiSFVecMass*complex<double> (1.,0.);
 

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

void  D0ToPiPiSKDec::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){

  if (_piPiSASHyp){
    std::map<std::string, double>::const_iterator itbFac;
    for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){ 
      fitVal.otherParams[itbFac->first]=itbFac->second;
      fitErr.otherParams[itbFac->first]=1.0;
    }
    
    std::map<std::string, double>::const_iterator itfFac;
    for(itfFac=_currentfProdFactorMap.begin();itfFac!=_currentfProdFactorMap.end(); ++itfFac){ 
      fitVal.otherParams[itfFac->first]=itfFac->second;
      fitErr.otherParams[itfFac->first]=1.0;
    }
    
    fitVal.otherParams[_name+"S0prodPosNeg"]=_currentS0Val; 
    fitErr.otherParams[_name+"S0prodPosNeg"]=1.0;
  }

}

void D0ToPiPiSKDec::print(std::ostream& os) const{
  return; //dummy
}

void D0ToPiPiSKDec::initialize(){
  if(_name==_ASKey){
    _piPiSASHyp=true;

    _currentbFactorMap[_name+"b_pole1Mag"]=1.;
    _currentbFactorMap[_name+"b_pole1Phi"]=0.;
    _currentbFactorMap[_name+"b_pole2Mag"]=1.;
    _currentbFactorMap[_name+"b_pole2Phi"]=0.;
    _currentbFactorMap[_name+"b_pole3Mag"]=1.;
    _currentbFactorMap[_name+"b_pole3Phi"]=0.;
    _currentbFactorMap[_name+"b_pole4Mag"]=1.;
    _currentbFactorMap[_name+"b_pole4Phi"]=0.;
    _currentbFactorMap[_name+"b_pole5Mag"]=1.;
    _currentbFactorMap[_name+"b_pole5Phi"]=0.;
    
    
    _currentfProdFactorMap[_name+"fprod_pipiMag"]=1.;
    _currentfProdFactorMap[_name+"fprod_pipiPhi"]=0.;
    _currentfProdFactorMap[_name+"fprod_KKMag"]=1.;
    _currentfProdFactorMap[_name+"fprod_KKPhi"]=0.;
    _currentfProdFactorMap[_name+"fprod_4piMag"]=1.;
    _currentfProdFactorMap[_name+"fprod_4piPhi"]=0.;
    _currentfProdFactorMap[_name+"fprod_etaetaMag"]=1.;
    _currentfProdFactorMap[_name+"fprod_etaetaPhi"]=0.;
    _currentfProdFactorMap[_name+"fprod_etaetapMag"]=1.;
    _currentfProdFactorMap[_name+"fprod_etaetapPhi"]=0.;

      complex<double> fProdPiPi= _currentfProdFactorMap[_name+"fprod_pipiMag"]*(complex<double> (cos(_currentbFactorMap[_name+"fprod_pipiPhi"]), sin(_currentbFactorMap[_name+"fprod_pipiPhi"])));
      complex<double> fProdKK= _currentfProdFactorMap[_name+"fprod_KKMag"]*(complex<double> (cos(_currentbFactorMap[_name+"fprod_KKPhi"]), sin(_currentbFactorMap[_name+"fprod_KKPhi"])));
      complex<double> fProd4Pi= _currentfProdFactorMap[_name+"fprod_4piMag"]*(complex<double> (cos(_currentbFactorMap[_name+"fprod_4piPhi"]), sin(_currentbFactorMap[_name+"fprod_4piPhi"])));
      complex<double> fProdEtaEta= _currentfProdFactorMap[_name+"fprod_etaetaMag"]*(complex<double> (cos(_currentbFactorMap[_name+"fprod_etaetaPhi"]), sin(_currentbFactorMap[_name+"fprod_etaetaPhi"])));
      complex<double> fProdEtaEtap= _currentfProdFactorMap[_name+"fprod_etaetapMag"]*(complex<double> (cos(_currentbFactorMap[_name+"fprod_etaetapPhi"]), sin(_currentbFactorMap[_name+"fprod_etaetapPhi"])));

      _pipiSFVec->updateFprod (0, fProdPiPi); 
      _pipiSFVec->updateFprod (1, fProdKK);
      _pipiSFVec->updateFprod (2, fProd4Pi);
      _pipiSFVec->updateFprod (3, fProdEtaEta);
      _pipiSFVec->updateFprod (4, fProdEtaEtap);
  }
  else{
    Alert << "PiPiS wave doesn't exist for key\t" << _name <<endmsg;
    exit(0);
  }
  Info << "hypothesis\t" << _name << "\t found" << endmsg;
}

bool D0ToPiPiSKDec::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=false;

  if (_piPiSASHyp){
    std::map<std::string, double>::const_iterator itbFac;
    for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){
      double currentbFactor=theParamVal.otherParams[itbFac->first];
      if ( fabs(currentbFactor-itbFac->second) > 1.e-10){
	DebugMsg << "bFactor " << itbFac->first << ":\t" << "current: " << itbFac->second << "\tnew: " << currentbFactor << endmsg;
	_recalculate=true;
      } 
    }

    std::map<std::string, double>::const_iterator itfProdFac;
    for(itfProdFac=_currentfProdFactorMap.begin();itfProdFac!=_currentfProdFactorMap.end(); ++itfProdFac){
      double currentfFactor=theParamVal.otherParams[itfProdFac->first];
      if ( fabs(currentfFactor-itfProdFac->second) > 1.e-10){
	DebugMsg << "fProdFactor " << itfProdFac->first << ":\t" << "current: " << itfProdFac->second << "\tnew: " << currentfFactor << endmsg;
	_recalculate=true;
      } 
    }
    
    if( fabs(_currentS0Val-theParamVal.otherParams[_name+"S0prodPosNeg"]) >1.e-10){
      DebugMsg << "S0ProdFactor " << _name << "S0prod:\t" << "current: " << _currentS0Val << "\tnew: " << theParamVal.otherParams[_name+"S0prodPosNeg"] << endmsg;
      _recalculate=true;
    }    
  }
  
  if (_recalculate) Info << "Recalculate amplitude:\t" << _name << endmsg;
  
  return _recalculate;
}


void D0ToPiPiSKDec::updateFitParams(fitParamsNew& theParamVal){
 
 if (_piPiSASHyp){
   std::map<std::string, double>::iterator itbFac;
   
   for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){
     double currentbFactor=theParamVal.otherParams[itbFac->first];
     itbFac->second = currentbFactor;
   }

   //update _pipiSFVec
   complex<double> b_pole1=_currentbFactorMap[_name+"b_pole1Mag"]*complex<double>(cos(_currentbFactorMap[_name+"b_pole1Phi"]), sin(_currentbFactorMap[_name+"b_pole1Phi"]));
   complex<double> b_pole2=_currentbFactorMap[_name+"b_pole2Mag"]*complex<double>(cos(_currentbFactorMap[_name+"b_pole2Phi"]), sin(_currentbFactorMap[_name+"b_pole2Phi"])); 
   complex<double> b_pole3=_currentbFactorMap[_name+"b_pole3Mag"]*complex<double>(cos(_currentbFactorMap[_name+"b_pole3Phi"]), sin(_currentbFactorMap[_name+"b_pole3Phi"]));
   complex<double> b_pole4=_currentbFactorMap[_name+"b_pole4Mag"]*complex<double>(cos(_currentbFactorMap[_name+"b_pole4Phi"]), sin(_currentbFactorMap[_name+"b_pole4Phi"]));
   complex<double> b_pole5=_currentbFactorMap[_name+"b_pole5Mag"]*complex<double>(cos(_currentbFactorMap[_name+"b_pole5Phi"]), sin(_currentbFactorMap[_name+"b_pole5Phi"]));
   
   _pipiSFVec->updateBeta(0, b_pole1);
   _pipiSFVec->updateBeta(1, b_pole2);
   _pipiSFVec->updateBeta(2, b_pole3);
   _pipiSFVec->updateBeta(3, b_pole4);
   _pipiSFVec->updateBeta(4, b_pole5);

   std::map<std::string, double>::iterator itfFac;
   
   for(itfFac=_currentfProdFactorMap.begin();itfFac!=_currentfProdFactorMap.end(); ++itfFac){
     double currentfFactor=theParamVal.otherParams[itfFac->first];
     itfFac->second = currentfFactor;
   }

   complex<double> fProdPiPi=_currentfProdFactorMap[_name+"fprod_pipiMag"]*complex<double>(cos(_currentfProdFactorMap[_name+"fprod_pipiPhi"]), sin(_currentfProdFactorMap[_name+"fprod_pipiPhi"]));
   complex<double> fProdKK=_currentfProdFactorMap[_name+"fprod_KKMag"]*complex<double>(cos(_currentfProdFactorMap[_name+"fprod_KKPhi"]), sin(_currentfProdFactorMap[_name+"fprod_KKPhi"]));
   complex<double> fProd4Pi=_currentfProdFactorMap[_name+"fprod_4piMag"]*complex<double>(cos(_currentfProdFactorMap[_name+"fprod_4piPhi"]), sin(_currentfProdFactorMap[_name+"fprod_4piPhi"]));
   complex<double> fProdEtaEta=_currentfProdFactorMap[_name+"fprod_etaetaMag"]*complex<double>(cos(_currentfProdFactorMap[_name+"fprod_etaetaPhi"]), sin(_currentfProdFactorMap[_name+"fprod_etaetaPhi"]));
   complex<double> fProdEtaEtap=_currentfProdFactorMap[_name+"fprod_etaetapMag"]*complex<double>(cos(_currentfProdFactorMap[_name+"fprod_etaetapPhi"]), sin(_currentfProdFactorMap[_name+"fprod_etaetapPhi"]));

   _pipiSFVec->updateFprod(0, fProdPiPi);
   _pipiSFVec->updateFprod(1, fProdKK);
   _pipiSFVec->updateFprod(2, fProd4Pi);
   _pipiSFVec->updateFprod(3, fProdEtaEta);
   _pipiSFVec->updateFprod(4, fProdEtaEtap);

   _currentS0Val=theParamVal.otherParams[_name+"S0prodPosNeg"];
   _pipiSFVec->updateS0prod(_currentS0Val);
 
 }


}
