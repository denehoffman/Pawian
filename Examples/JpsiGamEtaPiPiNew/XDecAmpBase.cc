#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
// #include "PwaUtils/EvtDataBaseListNew.hh"
#include "Examples/JpsiGamEtaPiPiNew/JpsiGamEtaPiPiEventListNew.hh"

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
{
  initialize();
}

XDecAmpBase::~XDecAmpBase()
{
}

complex<double> XDecAmpBase::XdecAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal){
  complex<double> result(0.,0.);

  if (_piPiEtaHyp){
    result+=XToPiPiEtaAmp(lamX, theData, theParamVal.Mags[_piPiEtaKey], theParamVal.Phis[_piPiEtaKey]);
  }
  if(_a980piHyp){
    double a0_980Mass=theParamVal.Masses["a0_980"];
    double a0_980gPiEta=theParamVal.gFactors["a0_980gPiEta"];
    double a0_980gKK=theParamVal.gFactors["a0_980gKK"];
    result+=XToAPiFlatteAmp(lamX, theData, theParamVal.Mags[_a980piKey], theParamVal.Phis[_a980piKey], a0_980Mass, a0_980gPiEta, a0_980gKK);
  }
  if(_a2_1320piHyp){
    double a2_1320Mass=theParamVal.Masses["a2_1320"];
    double a2_1320Width=theParamVal.Widths["a2_1320"];
    result+=XToAPiBWAmp(lamX, 2, theData, theParamVal.Mags[_a2_1320piKey], theParamVal.Phis[_a2_1320piKey] , a2_1320Mass, a2_1320Width);
  }
  if(_sigmaEtaHyp){
    double sigmaMass=theParamVal.Masses["sigma"];
    double sigmaWidth=theParamVal.Widths["sigma"];
    result+=XToEtaFAmp(lamX, 0, theData, theParamVal.Mags[_sigmaEtaKey], theParamVal.Phis[_sigmaEtaKey],sigmaMass, sigmaWidth);
  }
  if(_f980etaHyp){
    double f0_980Mass=theParamVal.Masses["f0_980"];
    double f0_980gPiPi=theParamVal.gFactors["f0_980gPiPi"];
    double f0_980gKK=theParamVal.gFactors["f0_980gKK"];
    result+=XToFEtaFlatteAmp(lamX, theData, theParamVal.Mags[_f980etaKey], theParamVal.Phis[_f980etaKey], f0_980Mass, f0_980gPiPi, f0_980gKK);
  }
  if(_f2_1270etaHyp){
    double f2_1270Mass=theParamVal.Masses["f2_1270"];
    double f2_1270Width=theParamVal.Widths["f2_1270"];
    result+=XToEtaFAmp(lamX, 2, theData, theParamVal.Mags[_f2_1270etaKey], theParamVal.Phis[_f2_1270etaKey], f2_1270Mass, f2_1270Width);
  }

  complex<double> dynModel(1.,0.);

  if (!_massIndependent){
    double xMass=theParamVal.Masses[_name];
    double xWidth=theParamVal.Widths[_name];
    Vector4<double> p4EtaPiPi = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPipPim_HeliPsi];
    dynModel=BreitWigner(p4EtaPiPi, xMass, xWidth);
  }

  result *=dynModel;
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
					     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& etaToA980PiMag,
					     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& etaToA980PiPhi,
					     double a0_980Mass, double a0_980gPiEta,double a0_980gKK){
  complex<double> result(0.,0.);
  Vector4<double > p4EtaPiplus = theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi];
  Vector4<double > p4EtaPiminus= theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi];

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

  for ( itXMag=etaToA980PiMag.begin(); itXMag!=etaToA980PiMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=etaToA980PiPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
        
    complex<double> amp = theXMag*expiphiX*sqrt(2.*XState->L+1.)*
      (  conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAplusPiminus][_J_X][lamX][0])*Flatte(p4EtaPiplus , _decPairPiEta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK)+
	 conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAminusPiplus][_J_X][lamX][0])*Flatte(p4EtaPiminus, _decPairPiEta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK)
	 );
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
      conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XTofEta][_J_X][lamX][0])*
      Flatte(p4PiPi, _decPairPiPi, _decPairKK, f0_980Mass, f0_980gPiPi, f0_980gKK);
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
      amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)
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

complex<double> XDecAmpBase::XToAPiBWAmp(Spin lamX, Spin jA, EvtDataNew* theData, 
			      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& etaToAPiMag, 
			      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& etaToAPiPhi, 
			      double aMass, double aWidth){

  complex<double> result(0.,0.);

  Vector4<double > p4EtaPiplus=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPip_HeliPsi];
  Vector4<double > p4EtaPiminus=theData->FourVecsDec[enumJpsiGamEtaPiPi4V::EtaPim_HeliPsi];  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

  for ( itXMag=etaToAPiMag.begin(); itXMag!=etaToAPiMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=etaToAPiPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
    complex<double> amp(0.,0.);     
    
    for(Spin lamA = -jA; lamA <= jA; lamA++){
      if(fabs(lamA)> _J_X || fabs(lamA)>XState->S) continue;
      
      amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)*sqrt(2.*jA + 1.)
	*Clebsch(XState->L, 0, XState->S, lamA, _J_X, lamA)
	*Clebsch(jA, lamA, 0, 0, XState->S, lamA)
	*(  conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAplusPiminus][_J_X][lamX][lamA])
	    *BreitWignerBlattW(p4EtaPiplus, _massPi, _massEta, aMass, aWidth, jA)
	    *conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::Aplus][jA][lamA][0])+
	    conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::XToAminusPiplus][_J_X][lamX][lamA])
	    *BreitWignerBlattW(p4EtaPiminus, _massPi, _massEta, aMass, aWidth, jA)
	    *conj(theData->WignerDsDec[enumJpsiGamEtaPiPiDfunc::Aminus][jA][lamA][0])
	 );
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

    fitVal.Masses["sigma"]=0.5;
    fitErr.Masses["sigma"]=0.09;
    fitVal.Widths["sigma"]=0.3;
    fitErr.Widths["sigma"]=0.059;

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
    }

    if (it->compare(0, _a980piKey.size(), _a980piKey)== 0){
      Info << "hypothesis\t" << _a980piKey << "\t enabled" << endmsg;
      _a980piHyp=true;
    }

    else if (it->compare(0, _sigmaEtaKey.size(), _sigmaEtaKey)== 0){
      Info << "hypothesis\t" << _sigmaEtaKey << "\t enabled" << endmsg;
      _sigmaEtaHyp=true;
    }

    else if (it->compare(0, _f980etaKey.size(), _f980etaKey)== 0){
      Info << "hypothesis\t" << _f980etaKey << "\t enabled" << endmsg;
      _f980etaHyp=true;
    }

    else if (it->compare(0, _a2_1320piKey.size(), _a2_1320piKey)== 0){
      Info << "hypothesis\t" << _a2_1320piKey << "\t enabled" << endmsg;
      _a2_1320piHyp=true;
    }

    else if (it->compare(0, _f2_1270etaKey.size(), _f2_1270etaKey)== 0){
      Info << "hypothesis\t" << _f2_1270etaKey << "\t enabled" << endmsg;
      _f2_1270etaHyp=true;
    }

    else if (it->compare(0, _xBWKey.size(), _xBWKey) ==0){
      _massIndependent=false;
    }
    
  }

}
