#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
// #include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiEventListNew.hh"
#include "PwaDynamics/FVectorPiPiS.hh"

XDecAmpBase::XDecAmpBase(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates, Spin spinX) :
  AbsXdecAmp(name, hypVec, spinX)
  ,_piPiEtaKey(name+"ToPiPiEta")
  ,_piPiEtaHyp(false)
  ,_a980piKey(name+"Toa980Pi")
  ,_a980piHyp(false)
  ,_sigmaEtaKey(name+"ToSigmaEta")
  ,_sigmaEtaHyp(false)
  ,_f980etaKey(name+"Tof980Eta")
  ,_f980etaHyp(false)
  ,_a2_1320piKey(name+"Toa2_1320Pi")
  ,_a2_1320piHyp(false)
  ,_f2_1270etaKey(name+"Tof2_1270Eta")
  ,_f2_1270etaHyp(false)
  ,_pipiSetaKey(name+"TopipiSEta")
  ,_pipiSetaHyp(false)
  ,_xBWKey(name+"BreitWigner")
  ,_massIndependent(true)
  ,_massPi(0.13957018)
  ,_massPi0(0.1349766)
  ,_massKplus(0.493677)
  ,_massK0(0.497614)
  ,_massEta(0.547853)
  ,_decPairKK(make_pair(_massKplus,_massK0))
  ,_decPairPiEta(make_pair(_massPi, _massEta))
  ,_decPairPi0Eta(make_pair(_massPi0, _massEta))
  ,_decPairPiPi(make_pair(_massPi, _massPi))
  ,_decPairPi0Pi0(make_pair(_massPi0, _massPi0))
  ,_theStatesPtr(theStates)
  ,_recalculatef2_1270(true)
  ,_recalculatea2_1320(true)
  ,_pipiSFVec(new FVectorPiPiS())
{
  initialize();
}

XDecAmpBase::~XDecAmpBase()
{
}

complex<double> XDecAmpBase::XdecAmp(Spin lamX, EvtDataNew* theData){
  int evtNo=theData->evtNo;
  if ( _cacheAmps && !_recalculate) return _cachedAmpMap[theData->evtNo][lamX];

  complex<double> result(0.,0.);

  if (_piPiEtaHyp){
    result+=XToPiPiEtaAmp(lamX, theData, _currentParamMagMap[_piPiEtaKey], _currentParamPhiMap[_piPiEtaKey]);
  }
  if (_pipiSetaHyp){
    result+=XToPiPiSEtaAmp(lamX, theData);
  }

  if(_a980piHyp){
    result+=XToAPiFlatteAmp(lamX, theData, _currentParamMagMap[_a980piKey], _currentParamPhiMap[_a980piKey], _currentMassMap["a0_980"], _currentgFactorMap["a0_980gPiEta"], _currentgFactorMap["a0_980gKK"]);
  }

  if(_a2_1320piHyp){
    if (_cacheAmps && !_recalculatea2_1320) result+=_cachedAmpa2_1320[evtNo][lamX];
    else {
      complex<double> tmpResult=XToAPiBWAmp(lamX, 2, theData, _currentParamMagMap[_a2_1320piKey], _currentParamPhiMap[_a2_1320piKey] , _currentMassMap["a2_1320"], _currentWidthMap["a2_1320"]);
      result+=tmpResult;
      if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	  _cachedAmpa2_1320[evtNo][lamX]=tmpResult;
#ifdef _OPENMP
	}
#endif
      }
    }
  }
  
  if(_sigmaEtaHyp){
    result+=XToEtaFAmp(lamX, 0, theData, _currentParamMagMap[_sigmaEtaKey], _currentParamPhiMap[_sigmaEtaKey],_currentMassMap["sigma"], _currentWidthMap["sigma"]);
  }
  if(_f980etaHyp){
    result+=XToFEtaFlatteAmp(lamX, theData, _currentParamMagMap[_f980etaKey], _currentParamPhiMap[_f980etaKey], _currentMassMap["f0_980"], _currentgFactorMap["f0_980gPiPi"], _currentgFactorMap["f0_980gKK"]);
  }
  if(_f2_1270etaHyp){
    if (_cacheAmps && !_recalculatef2_1270) result+=_cachedAmpf2_1270[evtNo][lamX];
    else {
      complex<double> tmpResult=XToEtaFAmp(lamX, 2, theData, _currentParamMagMap[_f2_1270etaKey], _currentParamPhiMap[_f2_1270etaKey], _currentMassMap["f2_1270"], _currentWidthMap["f2_1270"]);
      result+=tmpResult;
      if ( _cacheAmps){
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	  _cachedAmpf2_1270[evtNo][lamX]=tmpResult;
#ifdef _OPENMP
	}
#endif
      }
      
    }
  }

  complex<double> dynModel(1.,0.);

  if (!_massIndependent){
    Vector4<double> p4EtaPiPi = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPipPim_HeliPsi];
    dynModel=BreitWigner(p4EtaPiPi, _currentXMass, _currentXWidth);
  }

  result *=dynModel;

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

complex<double> XDecAmpBase::XToPiPiEtaAmp(Spin lamX, EvtDataNew* theData, 
					   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaMag,
					   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaPhi){

  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

   for ( itXMag=XToPiPiEtaMag.begin(); itXMag!=XToPiPiEtaMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToPiPiEtaPhi[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));

     result+= theXMag*expiphiX*sqrt(2.*XState->L+1.)
       *Clebsch(XState->L, 0, XState->S, 0, XState->J, 0)
       *Clebsch(0, 0, 0, 0, XState->S, 0);

   }
   result*=conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::etapipi][_J_X][lamX][0]);
   return result;
}



complex<double> XDecAmpBase::XToAPiFlatteAmp(Spin lamX, EvtDataNew* theData, 
					     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToA980PiMag,
					     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToA980PiPhi,
					     double a0_980Mass, double a0_980gPiEta,double a0_980gKK){
  complex<double> result(0.,0.);
  Vector4<double > p4EtaPiplus = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi];
  Vector4<double > p4EtaPiminus= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi];

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

  for ( itXMag=XToA980PiMag.begin(); itXMag!=XToA980PiMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=XToA980PiPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));

    complex<double> amp = theXMag*expiphiX*sqrt(2.*XState->L+1.)*
      (  conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAplusPiminus][_J_X][lamX][0])
	 //*sqrt(theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi].M()/theData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPip].P())
	 *Flatte(p4EtaPiplus , _decPairPiEta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK)
	 +
	 conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAminusPiplus][_J_X][lamX][0])
	 //*sqrt(theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi].M()/theData->FourVecsDec[enumJpsiGamEtaPiPi4V::Eta_HeliEtaPim].P())
	 *Flatte(p4EtaPiminus, _decPairPiEta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK)
	 )*0.5;
    result+= amp;
  }
  
  return result;
}

complex<double> XDecAmpBase::XToFEtaFlatteAmp(Spin lamX, EvtDataNew* theData, 
					      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XTof980etaMag,
					      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XTof980etaPhi,
					      double f0_980Mass, double f0_980gPiPi,double f0_980gKK){
  complex<double> result(0.,0.);
  Vector4<double > p4PiPi=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliPsi]; 
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

  for ( itXMag=XTof980etaMag.begin(); itXMag!=XTof980etaMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=XTof980etaPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
     
    complex<double> amp = theXMag*expiphiX*sqrt(2.*XState->L+1.)*
      conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XTofEta][_J_X][lamX][0])
      *Flatte(p4PiPi, _decPairPiPi, _decPairKK, f0_980Mass, f0_980gPiPi, f0_980gKK);
    result+= amp;
  }

  return result;
}


complex<double> XDecAmpBase::XToEtaFAmp(Spin lamX, Spin jf, EvtDataNew* theData, 
					std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFMag, 
					std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFPhi,
					double fMass, double fWidth){
  
  complex<double> result(0.,0.);
  
  Vector4<double > p4PiPi=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliPsi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;
  
  for ( itXMag=XToEtaFMag.begin(); itXMag!=XToEtaFMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=XToEtaFPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
    complex<double> amp(0.,0.);     
    for(Spin lamf = -jf; lamf <= jf; lamf++){
      if( fabs(lamf)> _J_X || fabs(lamf)>XState->S) continue;
      amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)*sqrt(2.*jf+1.)
	*Clebsch(XState->L, 0, XState->S, lamf, _J_X, lamf)
	*Clebsch(jf, lamf, 0, 0, XState->S, lamf)
	*conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XTofEta][_J_X][lamX][lamf])
	*conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::fToPiPi][jf][lamf][0])
	*BreitWignerBlattW(p4PiPi, _massPi, _massPi, fMass, fWidth, jf);
    }
    result+= amp;
  }
  return result;
}


complex<double> XDecAmpBase::XToPiPiSEtaAmp(Spin lamX, EvtDataNew* theData){
  complex<double> result(0.,0.);
  Vector4<double > p4PiPi=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::PipPim_HeliPsi];

  complex<double> pipiSFVecMass;
#ifdef _OPENMP
#pragma omp critical
  {
#endif
    _pipiSFVec->evalMatrix(p4PiPi.M());
    pipiSFVecMass=(*_pipiSFVec)[0];
#ifdef _OPENMP
  }
#endif  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > pipiSEtaMag=_currentParamMagMap[_pipiSetaKey];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > pipiSEtaPhi=_currentParamPhiMap[_pipiSetaKey];

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;  

  for ( itXMag=pipiSEtaMag.begin(); itXMag!=pipiSEtaMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=pipiSEtaPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
     
    complex<double> amp = theXMag*expiphiX*sqrt(2.*XState->L+1.)*
      conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XTofEta][_J_X][lamX][0])*pipiSFVecMass;
    result+= amp;
  }

  return result;
}



complex<double> XDecAmpBase::XToAPiBWAmp(Spin lamX, Spin jA, EvtDataNew* theData, 
			      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiMag, 
			      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiPhi, 
			      double aMass, double aWidth){

  complex<double> result(0.,0.);

  Vector4<double > p4EtaPiplus=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi];
  Vector4<double > p4EtaPiminus=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi];  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

  for ( itXMag=XToAPiMag.begin(); itXMag!=XToAPiMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=XToAPiPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
    complex<double> amp(0.,0.);     
    
    for(Spin lamA = -jA; lamA <= jA; lamA++){
      if(fabs(lamA)> _J_X || fabs(lamA)>XState->S){ 
	continue;
      }
      amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)*sqrt(2.*jA+1.)
	*Clebsch(XState->L, 0, XState->S, lamA, _J_X, lamA)
	*Clebsch(jA, lamA, 0, 0, XState->S, lamA)
	*(  conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAplusPiminus][_J_X][lamX][lamA])
	    *BreitWignerBlattW(p4EtaPiplus, _massPi, _massEta, aMass, aWidth, jA)
	    *conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::Aplus][jA][lamA][0])+
	    conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAminusPiplus][_J_X][lamX][lamA])
	    *BreitWignerBlattW(p4EtaPiminus, _massPi, _massEta, aMass, aWidth, jA)
	    *conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::Aminus][jA][lamA][0])
	 )*0.5;
    }
    result+= amp;
  }
  return result;
}


void  XDecAmpBase::getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr){
  
  if (_piPiEtaHyp){
    std::vector< boost::shared_ptr<const JPCLS> > PiPiEtaStates;
    if(_J_X==0) PiPiEtaStates=_theStatesPtr->EtaToa0PiStates();
    else if(_J_X==1) PiPiEtaStates=_theStatesPtr->F1Toa0PiStates();
    else if(_J_X==2) PiPiEtaStates=_theStatesPtr->Eta2Toa0PiStates();  
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;
    
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    for(itLS=PiPiEtaStates.begin(); itLS!=PiPiEtaStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_piPiEtaKey]=currentMagValMap;
    fitVal.Phis[_piPiEtaKey]=currentPhiValMap;
    fitErr.Mags[_piPiEtaKey]=currentMagErrMap;
    fitErr.Phis[_piPiEtaKey]=currentPhiErrMap;
  }

  if(_a980piHyp){
    std::vector< boost::shared_ptr<const JPCLS> > a0PiStates;
    if(_J_X==0) a0PiStates=_theStatesPtr->EtaToa0PiStates();
    else if(_J_X==1) a0PiStates=_theStatesPtr->F1Toa0PiStates();
    else if(_J_X==2) a0PiStates=_theStatesPtr->Eta2Toa0PiStates();  

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=a0PiStates.begin(); itLS!=a0PiStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_a980piKey]=currentMagValMap;
    fitVal.Phis[_a980piKey]=currentPhiValMap;
    fitErr.Mags[_a980piKey]=currentMagErrMap;
    fitErr.Phis[_a980piKey]=currentPhiErrMap;

    fitVal.Masses["a0_980"]=0.98;
    fitErr.Masses["a0_980"]=0.03;
    fitVal.gFactors["a0_980gPiEta"]=1.;
    fitErr.gFactors["a0_980gPiEta"]=0.19;
    fitVal.gFactors["a0_980gKK"]=1.;
    fitErr.gFactors["a0_980gKK"]=0.19;
  }

  if(_sigmaEtaHyp){
    std::vector< boost::shared_ptr<const JPCLS> > f0EtaStates;
    if(_J_X==0) f0EtaStates=_theStatesPtr->EtaTof0EtaStates();
    else if(_J_X==1) f0EtaStates=_theStatesPtr->F1Tof0EtaStates();
    else if(_J_X==2) f0EtaStates=_theStatesPtr->Eta2Tof0EtaStates(); 

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=f0EtaStates.begin(); itLS!=f0EtaStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_sigmaEtaKey]=currentMagValMap;
    fitVal.Phis[_sigmaEtaKey]=currentPhiValMap;
    fitErr.Mags[_sigmaEtaKey]=currentMagErrMap;
    fitErr.Phis[_sigmaEtaKey]=currentPhiErrMap;

    fitVal.Masses["sigma"]=0.475;
    fitErr.Masses["sigma"]=0.015;
    fitVal.Widths["sigma"]=0.550;
    fitErr.Widths["sigma"]=0.03;

  }

  if(_f980etaHyp){
    std::vector< boost::shared_ptr<const JPCLS> > f0EtaStates;
    if(_J_X==0) f0EtaStates=_theStatesPtr->EtaTof0EtaStates();
    else if(_J_X==1) f0EtaStates=_theStatesPtr->F1Tof0EtaStates();
    else if(_J_X==2) f0EtaStates=_theStatesPtr->Eta2Tof0EtaStates(); 

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=f0EtaStates.begin(); itLS!=f0EtaStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_f980etaKey]=currentMagValMap;
    fitVal.Phis[_f980etaKey]=currentPhiValMap;
    fitErr.Mags[_f980etaKey]=currentMagErrMap;
    fitErr.Phis[_f980etaKey]=currentPhiErrMap;

    fitVal.Masses["f0_980"]=0.98;
    fitErr.Masses["f0_980"]=0.03;
    fitVal.gFactors["f0_980gPiPi"]=0.7;
    fitErr.gFactors["f0_980gPiPi"]=0.2;
    fitVal.gFactors["f0_980gKK"]=1.2;
    fitErr.gFactors["f0_980gKK"]=0.2;
  }

  if(_a2_1320piHyp){
    std::vector< boost::shared_ptr<const JPCLS> > a2PiStates;
    if(_J_X==0) a2PiStates=_theStatesPtr->EtaToa2PiStates();
    else if(_J_X==1) a2PiStates=_theStatesPtr->F1Toa2PiStates();
    else if(_J_X==2) a2PiStates=_theStatesPtr->Eta2Toa2PiStates();  

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=a2PiStates.begin(); itLS!=a2PiStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_a2_1320piKey]=currentMagValMap;
    fitVal.Phis[_a2_1320piKey]=currentPhiValMap;
    fitErr.Mags[_a2_1320piKey]=currentMagErrMap;
    fitErr.Phis[_a2_1320piKey]=currentPhiErrMap;

    fitVal.Masses["a2_1320"]=1.32;
    fitErr.Masses["a2_1320"]=0.02;
    fitVal.Widths["a2_1320"]=0.1;
    fitErr.Widths["a2_1320"]=0.02;
  }

  if(_f2_1270etaHyp){
    std::vector< boost::shared_ptr<const JPCLS> > f2EtaStates;
    if(_J_X==0) f2EtaStates=_theStatesPtr->EtaTof2EtaStates();
    else if(_J_X==1) f2EtaStates=_theStatesPtr->F1Tof2EtaStates();
    else if(_J_X==2) f2EtaStates=_theStatesPtr->Eta2Tof2EtaStates(); 

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=f2EtaStates.begin(); itLS!=f2EtaStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }

    fitVal.Mags[_f2_1270etaKey]=currentMagValMap;
    fitVal.Phis[_f2_1270etaKey]=currentPhiValMap;
    fitErr.Mags[_f2_1270etaKey]=currentMagErrMap;
    fitErr.Phis[_f2_1270etaKey]=currentPhiErrMap;

    fitVal.Masses["f2_1270"]=1.27;
    fitErr.Masses["f2_1270"]=0.02;
    fitVal.Widths["f2_1270"]=0.185;
    fitErr.Widths["f2_1270"]=0.02;

  }

  if(_pipiSetaHyp){
    std::vector< boost::shared_ptr<const JPCLS> > pipiSEtaStates;
    if(_J_X==0) pipiSEtaStates=_theStatesPtr->EtaTof0EtaStates();
    else if(_J_X==1) pipiSEtaStates=_theStatesPtr->F1Tof0EtaStates();
    else if(_J_X==2) pipiSEtaStates=_theStatesPtr->Eta2Tof0EtaStates(); 

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiValMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentMagErrMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > currentPhiErrMap;

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itLS;
    
    for(itLS=pipiSEtaStates.begin(); itLS!=pipiSEtaStates.end(); ++itLS){
      currentMagValMap[*itLS]=0.2;
      currentPhiValMap[*itLS]=0.;
      currentMagErrMap[*itLS]=0.8;
      currentPhiErrMap[*itLS]=0.3;
    }
    fitVal.Mags[_pipiSetaKey]=currentMagValMap;
    fitVal.Phis[_pipiSetaKey]=currentPhiValMap;
    fitErr.Mags[_pipiSetaKey]=currentMagErrMap;
    fitErr.Phis[_pipiSetaKey]=currentPhiErrMap;

    std::map<std::string, double>::const_iterator itbFac;
    for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){ 
      fitVal.otherParams[itbFac->first]=itbFac->second;
    }

  }


  if (!_massIndependent){
    size_t pos=_name.find("_");
    std::string massMeVString=_name.substr(pos+1); 
    stringstream massMeVStrStream(massMeVString);
    int MassMeV;
    massMeVStrStream >> MassMeV;
    double MassGeV= ( (double) MassMeV)/1000.;

    fitVal.Masses[_name]=MassGeV;
    fitErr.Masses[_name]=0.01;
    fitVal.Widths[_name]=0.2;
    fitErr.Widths[_name]=0.02;
  }
}

void XDecAmpBase::print(std::ostream& os) const{
  return; //dummy
}

void XDecAmpBase::initialize(){
  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){

    if (it->compare(0, _piPiEtaKey.size(), _piPiEtaKey)== 0){
      Info << "hypothesis\t" << _piPiEtaKey << "\t enabled" << endmsg;
      _piPiEtaHyp=true;
      _enabledAmpKeys.push_back(_piPiEtaKey);
    }

    if (it->compare(0, _a980piKey.size(), _a980piKey)== 0){
      Info << "hypothesis\t" << _a980piKey << "\t enabled" << endmsg;
      _a980piHyp=true;
      _enabledAmpKeys.push_back(_a980piKey);
      _enabledMassKeys.push_back("a0_980");
      _enabledFactorKeys.push_back("a0_980gPiEta");
      _enabledFactorKeys.push_back("a0_980gKK");
    }

    else if (it->compare(0, _sigmaEtaKey.size(), _sigmaEtaKey)== 0){
      Info << "hypothesis\t" << _sigmaEtaKey << "\t enabled" << endmsg;
      _sigmaEtaHyp=true;
      _enabledAmpKeys.push_back(_sigmaEtaKey);
      _enabledMassKeys.push_back("sigma");
      _enabledWidthKeys.push_back("sigma");
    }

    else if (it->compare(0, _f980etaKey.size(), _f980etaKey)== 0){
      Info << "hypothesis\t" << _f980etaKey << "\t enabled" << endmsg;
      _f980etaHyp=true;
      _enabledAmpKeys.push_back(_f980etaKey);
      _enabledMassKeys.push_back("f0_980");
      _enabledFactorKeys.push_back("f0_980gPiPi");
      _enabledFactorKeys.push_back("f0_980gKK");
    }

    else if (it->compare(0, _a2_1320piKey.size(), _a2_1320piKey)== 0){
      Info << "hypothesis\t" << _a2_1320piKey << "\t enabled" << endmsg;
      _a2_1320piHyp=true;
      _enabledAmpKeys.push_back(_a2_1320piKey);
      _enabledMassKeys.push_back("a2_1320");
      _enabledWidthKeys.push_back("a2_1320");
    }

    else if (it->compare(0, _f2_1270etaKey.size(), _f2_1270etaKey)== 0){
      Info << "hypothesis\t" << _f2_1270etaKey << "\t enabled" << endmsg;
      _f2_1270etaHyp=true;
      _enabledAmpKeys.push_back(_f2_1270etaKey);
      _enabledMassKeys.push_back("f2_1270");
      _enabledWidthKeys.push_back("f2_1270");
    }
    
    else if (it->compare(0, _pipiSetaKey.size(), _pipiSetaKey)== 0){
      Info << "hypothesis\t" << _pipiSetaKey << "\t enabled" << endmsg;
      _pipiSetaHyp=true;
      _enabledAmpKeys.push_back(_pipiSetaKey);
      _currentbFactorMap[_pipiSetaKey+"b_pole1Mag"]=1.;
      _currentbFactorMap[_pipiSetaKey+"b_pole1Phi"]=0.;
      _currentbFactorMap[_pipiSetaKey+"b_pole2Mag"]=1.;
      _currentbFactorMap[_pipiSetaKey+"b_pole2Phi"]=0.;
      _currentbFactorMap[_pipiSetaKey+"b_pole3Mag"]=1.;
      _currentbFactorMap[_pipiSetaKey+"b_pole3Phi"]=0.;
      _currentbFactorMap[_pipiSetaKey+"b_pole4Mag"]=1.;
      _currentbFactorMap[_pipiSetaKey+"b_pole4Phi"]=0.;
      _currentbFactorMap[_pipiSetaKey+"b_pole5Mag"]=1.;
      _currentbFactorMap[_pipiSetaKey+"b_pole5Phi"]=0.;
       
      complex<double> woFprod=0.;
      for (int i=0; i<5; ++i) _pipiSFVec->updateFprod (i, woFprod);

    }
    else if (it->compare(0, _xBWKey.size(), _xBWKey) ==0){
      _massIndependent=false;
    }
    
  }

}

void XDecAmpBase::checkRecalculation(fitParamsNew& theParamVal){
  _recalculate=false;

  std::vector<std::string>::const_iterator itKeys;
  for ( itKeys=_enabledAmpKeys.begin(); itKeys!=_enabledAmpKeys.end(); ++itKeys){ 
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > theMags=theParamVal.Mags[(*itKeys)];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > thePhis=theParamVal.Phis[(*itKeys)];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for ( itMag=theMags.begin(); itMag!=theMags.end(); ++itMag){
      boost::shared_ptr<const JPCLS> XState=itMag->first;
      double theMag=itMag->second;
      double thePhi=thePhis[XState];
      if ( fabs(theMag - _currentParamMagMap[(*itKeys)][XState])  > 1.e-10 ){
        _recalculate=true;
        // _currentParamMags[XState]=theMag;
        DebugMsg << setprecision (8) << "Difference Mag " << XState->name() << ":\t" << "current: " << _currentParamMagMap[(*itKeys)][XState] << "\tnew: " << theMag << endmsg;
      }
      if ( fabs(thePhi - _currentParamPhiMap[(*itKeys)][XState])  > 1.e-10 ){
        _recalculate=true;
        DebugMsg  << setprecision (8) << "Difference Phi " << XState->name() << ":\t" << "current: " << _currentParamPhiMap[(*itKeys)][XState] << "\tnew: " << thePhi << endmsg;
      }
    }
  }
 
 for ( itKeys=_enabledMassKeys.begin(); itKeys!=_enabledMassKeys.end(); ++itKeys){
   double currentMass=theParamVal.Masses[(*itKeys)];
   if ( fabs(currentMass-_currentMassMap[(*itKeys)]) > 1.e-10){
     DebugMsg << "Mass " << (*itKeys) << ":\t" << "current: " << _currentMassMap[(*itKeys)] << "\tnew: " << currentMass << endmsg;
     _recalculate=true;
   }   
 }

 for ( itKeys=_enabledWidthKeys.begin(); itKeys!=_enabledWidthKeys.end(); ++itKeys){
   double currentWidth=theParamVal.Widths[(*itKeys)];
   if ( fabs(currentWidth-_currentWidthMap[(*itKeys)]) > 1.e-10){
     DebugMsg << "Width " << (*itKeys) << ":\t" << "current: " << _currentWidthMap[(*itKeys)] << "\tnew: " << currentWidth << endmsg;
     _recalculate=true;
   }   
 }

 for ( itKeys=_enabledFactorKeys.begin(); itKeys!=_enabledFactorKeys.end(); ++itKeys){
   double currentgFactor=theParamVal.gFactors[(*itKeys)];
   if ( fabs(currentgFactor-_currentgFactorMap[(*itKeys)]) > 1.e-10){
     DebugMsg << "gFactor " << (*itKeys) << ":\t" << "current: " << _currentgFactorMap[(*itKeys)] << "\tnew: " << currentgFactor << endmsg;
     _recalculate=true;
   }   
 } 

 if (_pipiSetaHyp){
   std::map<std::string, double>::const_iterator itbFac;

   for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){
     double currentbFactor=theParamVal.otherParams[itbFac->first];
     if ( fabs(currentbFactor-itbFac->second) > 1.e-10){
       DebugMsg << "bFactor " << itbFac->first << ":\t" << "current: " << itbFac->second << "\tnew: " << currentbFactor << endmsg;
       _recalculate=true;
     } 
   }

 }


  if (!_massIndependent){
    double xMass=theParamVal.Masses[_name];
    if ( fabs(xMass-_currentXMass) > 1.e-10){ 
      _recalculate=true;
    }

    double xWidth=theParamVal.Widths[_name];
    if ( fabs(xWidth - _currentXWidth) > 1.e-10){
      _recalculate=true;
    }
  }
  if (_recalculate) Info << "Recalculate amplitude:\t" << _name << endmsg;

  //a2_1320 amp
  _recalculatea2_1320=false;
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > theMags=theParamVal.Mags[_a2_1320piKey];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > thePhis=theParamVal.Phis[_a2_1320piKey];
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
  for ( itMag=theMags.begin(); itMag!=theMags.end(); ++itMag){
    boost::shared_ptr<const JPCLS> XState=itMag->first;
    double theMag=itMag->second;
    double thePhi=thePhis[XState];
    if ( fabs(theMag - _currentParamMagMap[_a2_1320piKey][XState])  > 1.e-10 ){
      _recalculatea2_1320=true;
      DebugMsg << setprecision (8) << "Difference Mag " << XState->name() << ":\t" << "current: " << _currentParamMagMap[_a2_1320piKey][XState] << "\tnew: " << theMag << endmsg;
    }
    if ( fabs(thePhi - _currentParamPhiMap[_a2_1320piKey][XState])  > 1.e-10 ){
      _recalculatea2_1320=true;
      DebugMsg  << setprecision (8) << "Difference Phi " << XState->name() << ":\t" << "current: " << _currentParamPhiMap[_a2_1320piKey][XState] << "\tnew: " << thePhi << endmsg;
    }
  }

  double currentMass=theParamVal.Masses["a2_1320"];
  if ( fabs(currentMass-_currentMassMap["a2_1320"]) > 1.e-10){
    DebugMsg << "Mass a2_1320:\t" << "current: " << _currentMassMap["a2_1320"] << "\tnew: " << currentMass << endmsg;
    _recalculatea2_1320=true;
  }

  double currentWidth=theParamVal.Widths["a2_1320"];
  if ( fabs(currentWidth-_currentWidthMap["a2_1320"]) > 1.e-10){
    DebugMsg << "Width a2_1320:\t" << "current: " << _currentWidthMap["a2_1320"] << "\tnew: " << currentWidth << endmsg;
    _recalculatea2_1320=true;
  }

  if (_recalculatea2_1320) DebugMsg << "Recalculate a2_1320 amplitude in:\t" << _name << endmsg;

  //f2_1270 amp
  _recalculatef2_1270=false;
  theMags=theParamVal.Mags[_f2_1270etaKey];
  thePhis=theParamVal.Phis[_f2_1270etaKey];
  
  for ( itMag=theMags.begin(); itMag!=theMags.end(); ++itMag){
    boost::shared_ptr<const JPCLS> XState=itMag->first;
    double theMag=itMag->second;
    double thePhi=thePhis[XState];
    if ( fabs(theMag - _currentParamMagMap[_f2_1270etaKey][XState])  > 1.e-10 ){
      _recalculatef2_1270=true;
      DebugMsg << setprecision (8) << "Difference Mag " << XState->name() << ":\t" << "current: " << _currentParamMagMap[_f2_1270etaKey][XState] << "\tnew: " << theMag << endmsg;
    }
    if ( fabs(thePhi - _currentParamPhiMap[_f2_1270etaKey][XState])  > 1.e-10 ){
      _recalculatef2_1270=true;
      DebugMsg  << setprecision (8) << "Difference Phi " << XState->name() << ":\t" << "current: " << _currentParamPhiMap[_f2_1270etaKey][XState] << "\tnew: " << thePhi << endmsg;
    }
  }

  currentMass=theParamVal.Masses["f2_1270"];
  if ( fabs(currentMass-_currentMassMap["f2_1270"]) > 1.e-10){
    DebugMsg << setprecision (8) << "Mass f2_1270:\t" << "current: " << _currentMassMap["f2_1270"] << "\tnew: " << currentMass << endmsg;
    _recalculatef2_1270=true;
  }

  currentWidth=theParamVal.Widths["f2_1270"];
  if ( fabs(currentWidth-_currentWidthMap["f2_1270"]) > 1.e-10){
    DebugMsg << setprecision (8) << "Width f2_1270:\t" << "current: " << _currentWidthMap["f2_1270"] << "\tnew: " << currentWidth << endmsg;
    _recalculatef2_1270=true;
  }

  if (_recalculatef2_1270) Info << "Recalculate f2_1270 amplitude in:\t" << _name << endmsg;
}

void XDecAmpBase::updateFitParams(fitParamsNew& theParamVal){
  std::vector<std::string>::const_iterator itKeys;
  for ( itKeys=_enabledAmpKeys.begin(); itKeys!=_enabledAmpKeys.end(); ++itKeys){ 
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > theMags=theParamVal.Mags[(*itKeys)];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > thePhis=theParamVal.Phis[(*itKeys)];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for ( itMag=theMags.begin(); itMag!=theMags.end(); ++itMag){
      boost::shared_ptr<const JPCLS> XState=itMag->first;
      double theMag=itMag->second;
      double thePhi=thePhis[XState];
      _currentParamMagMap[(*itKeys)][XState]=theMag;
      _currentParamPhiMap[(*itKeys)][XState]=thePhi;
    }
  }

 for ( itKeys=_enabledMassKeys.begin(); itKeys!=_enabledMassKeys.end(); ++itKeys){
   double currentMass=theParamVal.Masses[(*itKeys)];
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

 if (_pipiSetaHyp){
   std::map<std::string, double>::iterator itbFac;
   
   for(itbFac=_currentbFactorMap.begin();itbFac!=_currentbFactorMap.end(); ++itbFac){
     double currentbFactor=theParamVal.otherParams[itbFac->first];
     itbFac->second = currentbFactor;
   }

   //update _pipiSFVec
   complex<double> b_pole1=_currentbFactorMap[_pipiSetaKey+"b_pole1Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSetaKey+"b_pole1Phi"]), sin(_currentbFactorMap[_pipiSetaKey+"b_pole1Phi"]));
   complex<double> b_pole2=_currentbFactorMap[_pipiSetaKey+"b_pole2Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSetaKey+"b_pole2Phi"]), sin(_currentbFactorMap[_pipiSetaKey+"b_pole2Phi"])); 
   complex<double> b_pole3=_currentbFactorMap[_pipiSetaKey+"b_pole3Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSetaKey+"b_pole3Phi"]), sin(_currentbFactorMap[_pipiSetaKey+"b_pole3Phi"]));
   complex<double> b_pole4=_currentbFactorMap[_pipiSetaKey+"b_pole4Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSetaKey+"b_pole4Phi"]), sin(_currentbFactorMap[_pipiSetaKey+"b_pole4Phi"]));
   complex<double> b_pole5=_currentbFactorMap[_pipiSetaKey+"b_pole5Mag"]*complex<double>(cos(_currentbFactorMap[_pipiSetaKey+"b_pole5Phi"]), sin(_currentbFactorMap[_pipiSetaKey+"b_pole5Phi"]));
   
   _pipiSFVec->updateBeta(0, b_pole1);
   _pipiSFVec->updateBeta(1, b_pole2);
   _pipiSFVec->updateBeta(2, b_pole3);
   _pipiSFVec->updateBeta(3, b_pole4);
   _pipiSFVec->updateBeta(4, b_pole5); 
 }

 if (!_massIndependent){
   double xMass=theParamVal.Masses[_name];
   _currentXMass= xMass;
   
   double xWidth=theParamVal.Widths[_name];
   _currentXWidth=xWidth;
   
 }
}
