#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Examples/D0ToKsPipPim/D0ToKsPipPimLh.hh"
#include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimEventList.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaUtils/FitParamsBaseNew.hh"
#include "PwaUtils/AbsXdecAmp.hh"
#include "PwaDynamics/FVectorPiPiS.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


D0ToKsPipPimLh::D0ToKsPipPimLh(boost::shared_ptr<const EvtDataBaseListNew> theEvtList, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates) :
  AbsLhNew(theEvtList)
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false) 
  ,_pipiSHyp(false)
  ,_pipiSKey("PiPiS")
  ,_phasespaceKey("Phasespace")
  ,_recalculatepipiS(true)
  ,_pipiSFVec(new FVectorPiPiS())
  ,_currentS0Val(0.)
{
  initializeHypothesis();
}

D0ToKsPipPimLh::D0ToKsPipPimLh( boost::shared_ptr<AbsLhNew> theLhPtr, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates) :
  AbsLhNew(theLhPtr->getEventList())
  ,_hypVec(hypVec)
  ,_theStatesPtr(theStates)
  ,_usePhasespace(false) 
  ,_pipiSHyp(false)
  ,_pipiSKey("PiPiS_")
  ,_phasespaceKey("Phasespace")
  ,_recalculatepipiS(true)
  ,_pipiSFVec(new FVectorPiPiS())
  ,_currentS0Val(0.)
{
  initializeHypothesis();  
}

D0ToKsPipPimLh::~D0ToKsPipPimLh()
{;
}

double D0ToKsPipPimLh::calcEvtIntensity(EvtDataNew* theData, fitParamsNew& theParamVal){

  double result=0.;

  complex<double> tmpAmp(0.,0.);

  if (_pipiSHyp){
    tmpAmp += D0ToPiPiSAmp(theData);
  }

  result=norm(tmpAmp);

  if(_usePhasespace){
    result = result + theParamVal.otherParams[_phasespaceKey];  
  }

  return result;  

}

complex<double> D0ToKsPipPimLh::D0ToPiPiSAmp(EvtDataNew* theData){
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

 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > pipiSMag=_currentParamMagMap[_pipiSKey];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > pipiSPhi=_currentParamPhiMap[_pipiSKey];

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;  

  for ( itXMag=pipiSMag.begin(); itXMag!=pipiSMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=pipiSPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
     
    complex<double> amp = theXMag*expiphiX*pipiSFVecMass;
    result+= amp;
  }

  return result;

}

void D0ToKsPipPimLh::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){ 

  if(_usePhasespace){
    fitVal.otherParams[_phasespaceKey]=0.02;
    fitErr.otherParams[_phasespaceKey]=0.015;
  }

  if(_pipiSHyp){
    std::vector< boost::shared_ptr<const JPCLS> > pipiSStates=_theStatesPtr->D0Tof0KStates();
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=pipiSStates.begin(); itLS!=pipiSStates.end(); ++itLS){
      currentMagValMap[*itLS]=1.;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }

    fitVal.Mags[_pipiSKey]=currentMagValMap;
    fitVal.Phis[_pipiSKey]=currentPhiValMap;
    fitErr.Mags[_pipiSKey]=currentMagErrMap;
    fitErr.Phis[_pipiSKey]=currentPhiErrMap;

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

    fitVal.otherParams[_pipiSKey+"S0prodPosNeg"]=_currentS0Val; 
    fitErr.otherParams[_pipiSKey+"S0prodPosNeg"]=1.0;
  }
}

void D0ToKsPipPimLh::print(std::ostream& os) const{

}


void  D0ToKsPipPimLh::initializeHypothesis(){

  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){

    if (it->compare(0, _pipiSKey.size(), _pipiSKey)== 0){
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _pipiSHyp=true;

      _enabledAmpKeys.push_back(_pipiSKey);
      _currentbFactorMap[_pipiSKey+"b_pole1Mag"]=1.;
      _currentbFactorMap[_pipiSKey+"b_pole1Phi"]=0.;
      _currentbFactorMap[_pipiSKey+"b_pole2Mag"]=1.;
      _currentbFactorMap[_pipiSKey+"b_pole2Phi"]=0.;
      _currentbFactorMap[_pipiSKey+"b_pole3Mag"]=1.;
      _currentbFactorMap[_pipiSKey+"b_pole3Phi"]=0.;
      _currentbFactorMap[_pipiSKey+"b_pole4Mag"]=1.;
      _currentbFactorMap[_pipiSKey+"b_pole4Phi"]=0.;
      _currentbFactorMap[_pipiSKey+"b_pole5Mag"]=1.;
      _currentbFactorMap[_pipiSKey+"b_pole5Phi"]=0.;


      _currentfProdFactorMap[_pipiSKey+"fprod_pipiMag"]=1.;
      _currentfProdFactorMap[_pipiSKey+"fprod_pipiPhi"]=0.;
      _currentfProdFactorMap[_pipiSKey+"fprod_KKMag"]=1.;
      _currentfProdFactorMap[_pipiSKey+"fprod_KKPhi"]=0.;
      _currentfProdFactorMap[_pipiSKey+"fprod_4piMag"]=1.;
      _currentfProdFactorMap[_pipiSKey+"fprod_4piPhi"]=0.;
      _currentfProdFactorMap[_pipiSKey+"fprod_etaetaMag"]=1.;
      _currentfProdFactorMap[_pipiSKey+"fprod_etaetaPhi"]=0.;
      _currentfProdFactorMap[_pipiSKey+"fprod_etaetapMag"]=1.;
      _currentfProdFactorMap[_pipiSKey+"fprod_etaetapPhi"]=0.;

      complex<double> fProdPiPi= _currentfProdFactorMap[_pipiSKey+"fprod_pipiMag"]*(complex<double> (cos(_currentbFactorMap[_pipiSKey+"fprod_pipiPhi"]), sin(_currentbFactorMap[_pipiSKey+"fprod_pipiPhi"])));
      complex<double> fProdKK= _currentfProdFactorMap[_pipiSKey+"fprod_KKMag"]*(complex<double> (cos(_currentbFactorMap[_pipiSKey+"fprod_KKPhi"]), sin(_currentbFactorMap[_pipiSKey+"fprod_KKPhi"])));
      complex<double> fProd4Pi= _currentfProdFactorMap[_pipiSKey+"fprod_4piMag"]*(complex<double> (cos(_currentbFactorMap[_pipiSKey+"fprod_4piPhi"]), sin(_currentbFactorMap[_pipiSKey+"fprod_4piPhi"])));
      complex<double> fProdEtaEta= _currentfProdFactorMap[_pipiSKey+"fprod_etaetaMag"]*(complex<double> (cos(_currentbFactorMap[_pipiSKey+"fprod_etaetaPhi"]), sin(_currentbFactorMap[_pipiSKey+"fprod_etaetaPhi"])));
      complex<double> fProdEtaEtap= _currentfProdFactorMap[_pipiSKey+"fprod_etaetapMag"]*(complex<double> (cos(_currentbFactorMap[_pipiSKey+"fprod_etaetapPhi"]), sin(_currentbFactorMap[_pipiSKey+"fprod_etaetapPhi"])));

      _pipiSFVec->updateFprod (0, fProdPiPi); 
      _pipiSFVec->updateFprod (1, fProdKK);
      _pipiSFVec->updateFprod (2, fProd4Pi);
      _pipiSFVec->updateFprod (3, fProdEtaEta);
      _pipiSFVec->updateFprod (4, fProdEtaEtap);
    }
    else if (it->compare(0, _phasespaceKey.size(), _phasespaceKey)== 0){
      Info << "hypothesis\t" << (*it) << "\t found" << endmsg;
      _usePhasespace=true;
    }
   
  } 

}



void D0ToKsPipPimLh::checkParamVariation(fitParamsNew& theParamVal){

 if (_pipiSHyp){

   _recalculatepipiS=false;
   std::map<std::string, double>::const_iterator itbFac;

   for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){
     double currentbFactor=theParamVal.otherParams[itbFac->first];
     if ( fabs(currentbFactor-itbFac->second) > 1.e-10){
       DebugMsg << "bFactor " << itbFac->first << ":\t" << "current: " << itbFac->second << "\tnew: " << currentbFactor << endmsg;
       _recalculatepipiS=true;
     } 
   }

   std::map<std::string, double>::const_iterator itfProdFac;
   for(itfProdFac=_currentfProdFactorMap.begin();itfProdFac!=_currentfProdFactorMap.end(); ++itfProdFac){
     double currentfFactor=theParamVal.otherParams[itfProdFac->first];
     if ( fabs(currentfFactor-itfProdFac->second) > 1.e-10){
       DebugMsg << "fProdFactor " << itfProdFac->first << ":\t" << "current: " << itfProdFac->second << "\tnew: " << currentfFactor << endmsg;
       _recalculatepipiS=true;
     } 
   }

   if( fabs(_currentS0Val-theParamVal.otherParams[_pipiSKey+"S0prodPosNeg"]) >1.e-10){
     DebugMsg << "S0ProdFactor " << _pipiSKey << "S0prod:\t" << "current: " << _currentS0Val << "\tnew: " << theParamVal.otherParams[_pipiSKey+"S0prodPosNeg"] << endmsg;
       _recalculatepipiS=true;
   }


   //for chaching: check pipiS amp  
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > theMags=theParamVal.Mags[_pipiSKey];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > thePhis=theParamVal.Phis[_pipiSKey];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for ( itMag=theMags.begin(); itMag!=theMags.end(); ++itMag){
      boost::shared_ptr<const JPCLS> XState=itMag->first;
      double theMag=itMag->second;
      double thePhi=thePhis[XState];
      if ( fabs(theMag - _currentParamMagMap[_pipiSKey][XState])  > 1.e-10 ){
        _recalculatepipiS=true;
        DebugMsg << setprecision (8) << "Difference Mag " << XState->name() << ":\t" << "current: " << _currentParamMagMap[_pipiSKey][XState] << "\tnew: " << theMag << endmsg;
      }
      if ( fabs(thePhi - _currentParamPhiMap[_pipiSKey][XState])  > 1.e-10 ){
        _recalculatepipiS=true;
        DebugMsg  << setprecision (8) << "Difference Phi " << XState->name() << ":\t" << "current: " << _currentParamPhiMap[_pipiSKey][XState] << "\tnew: " << thePhi << endmsg;
      }
    }
    if (_recalculatepipiS) DebugMsg << "Recalculate pipiS amplitude in:\t" << endmsg;
 }


  return;
}

void D0ToKsPipPimLh::cacheTheAmps(){
  return;
}

void D0ToKsPipPimLh::updateFitParams(fitParamsNew& theParamVal){
  std::vector<std::string>::const_iterator itKeys;
  for ( itKeys=_enabledAmpKeys.begin(); itKeys!=_enabledAmpKeys.end(); ++itKeys){ 
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess
 > theMags=theParamVal.Mags[(*itKeys)];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess
 > thePhis=theParamVal.Phis[(*itKeys)];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess
 >::iterator itMag;
    for ( itMag=theMags.begin(); itMag!=theMags.end(); ++itMag){
      boost::shared_ptr<const JPCLS> XState=itMag->first;
      double theMag=itMag->second;
      double thePhi=thePhis[XState];
      _currentParamMagMap[(*itKeys)][XState]=theMag;
      _currentParamPhiMap[(*itKeys)][XState]=thePhi;
    }
  }

for ( itKeys=_enabledMassKeys.begin(); itKeys!=_enabledMassKeys.end(); ++itKeys){   double currentMass=theParamVal.Masses[(*itKeys)];
   _currentMassMap[(*itKeys)]=currentMass;      
 }

 for ( itKeys=_enabledWidthKeys.begin(); itKeys!=_enabledWidthKeys.end(); ++itKeys){
   double currentWidth=theParamVal.Widths[(*itKeys)];
   _currentWidthMap[(*itKeys)]=currentWidth;
 }

 for ( itKeys=_enabledFactorKeys.begin(); itKeys!=_enabledFactorKeys.end(); ++itKeys){
   double currentgFactor=theParamVal.gFactors[(*itKeys)];
   _currentgFactorMap[(*itKeys)]=currentgFactor;
 }

 if (_pipiSHyp){
   std::map<std::string, double>::iterator itbFac;
   
   for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){
     double currentbFactor=theParamVal.otherParams[itbFac->first];
     itbFac->second = currentbFactor;
   }

   //update _pipiSFVec
   complex<double> b_pole1=_currentbFactorMap[_pipiSKey+"b_pole1Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSKey+"b_pole1Phi"]), sin(_currentbFactorMap[_pipiSKey+"b_pole1Phi"]));
   complex<double> b_pole2=_currentbFactorMap[_pipiSKey+"b_pole2Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSKey+"b_pole2Phi"]), sin(_currentbFactorMap[_pipiSKey+"b_pole2Phi"])); 
   complex<double> b_pole3=_currentbFactorMap[_pipiSKey+"b_pole3Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSKey+"b_pole3Phi"]), sin(_currentbFactorMap[_pipiSKey+"b_pole3Phi"]));
   complex<double> b_pole4=_currentbFactorMap[_pipiSKey+"b_pole4Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSKey+"b_pole4Phi"]), sin(_currentbFactorMap[_pipiSKey+"b_pole4Phi"]));
   complex<double> b_pole5=_currentbFactorMap[_pipiSKey+"b_pole5Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSKey+"b_pole5Phi"]), sin(_currentbFactorMap[_pipiSKey+"b_pole5Phi"]));
   
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

   complex<double> fProdPiPi=_currentfProdFactorMap[_pipiSKey+"fprod_pipiMag"]*complex<double>(cos(_currentfProdFactorMap[_pipiSKey+"fprod_pipiPhi"]), sin(_currentfProdFactorMap[_pipiSKey+"fprod_pipiPhi"]));
   complex<double> fProdKK=_currentfProdFactorMap[_pipiSKey+"fprod_KKMag"]*complex<double>(cos(_currentfProdFactorMap[_pipiSKey+"fprod_KKPhi"]), sin(_currentfProdFactorMap[_pipiSKey+"fprod_KKPhi"]));
   complex<double> fProd4Pi=_currentfProdFactorMap[_pipiSKey+"fprod_4piMag"]*complex<double>(cos(_currentfProdFactorMap[_pipiSKey+"fprod_4piPhi"]), sin(_currentfProdFactorMap[_pipiSKey+"fprod_4piPhi"]));
   complex<double> fProdEtaEta=_currentfProdFactorMap[_pipiSKey+"fprod_etaetaMag"]*complex<double>(cos(_currentfProdFactorMap[_pipiSKey+"fprod_etaetaPhi"]), sin(_currentfProdFactorMap[_pipiSKey+"fprod_etaetaPhi"]));
   complex<double> fProdEtaEtap=_currentfProdFactorMap[_pipiSKey+"fprod_etaetapMag"]*complex<double>(cos(_currentfProdFactorMap[_pipiSKey+"fprod_etaetapPhi"]), sin(_currentfProdFactorMap[_pipiSKey+"fprod_etaetapPhi"]));

   _pipiSFVec->updateFprod(0, fProdPiPi);
   _pipiSFVec->updateFprod(1, fProdKK);
   _pipiSFVec->updateFprod(2, fProd4Pi);
   _pipiSFVec->updateFprod(3, fProdEtaEta);
   _pipiSFVec->updateFprod(4, fProdEtaEtap);

   _currentS0Val=theParamVal.otherParams[_pipiSKey+"S0prodPosNeg"];
   _pipiSFVec->updateS0prod(_currentS0Val); 
 }
  return;
}
