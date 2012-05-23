#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiData.hh"

XDecAmpBase::XDecAmpBase(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates, Spin spinX) :
  AbsXdecAmp(name, hypVec)
  ,_a980piKey(name+"Toa980")
  ,_a980piHyp(false)
  ,_sigmaEtaKey(name+"ToSigmaEta")
  ,_sigmaEtaHyp(false)
  ,_massPi0(0.1349766)
  ,_massKplus(0.493677)
  ,_massK0(0.497614)
  ,_massEta(0.547853)
  ,_decPairKK(make_pair(_massKplus,_massK0))
  ,_decPairPi0Eta(make_pair(_massPi0, _massEta))
  ,_theStatesPtr(theStates)
  ,_spinX(spinX)
{
  initialize();
}

XDecAmpBase::~XDecAmpBase()
{
}



complex<double> XDecAmpBase::XToAPiFlatteAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal){
  complex<double> result(0.,0.);


  double a0_980Mass=theParamVal.Masses["a0_980"];
  double a0_980gPiEta=theParamVal.gFactors["a0_980gPiEta"];
  double a0_980gKK=theParamVal.gFactors["a0_980gKK"];

  Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
			       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
			       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
			       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());

  Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());


   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToA980PiMag=theParamVal.Mags[_a980piKey];

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToA980PiPhi=theParamVal.Phis[_a980piKey];


  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;


   for ( itXMag=etaToA980PiMag.begin(); itXMag!=etaToA980PiMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=etaToA980PiPhi[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));


     complex<double> amp = theXMag*expiphiX*sqrt(2.*XState->L+1.)*
       (  conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_XToAplusPiminusdec][_spinX][lamX][0])*Flatte(p4EtaPiplus, _decPairPi0Eta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK)+
	  conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_XToAminusPiplusdec][_spinX][lamX][0])*
	  Flatte(p4EtaPiminus, _decPairPi0Eta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK)
	 );
     result+= amp;
   }

  return result;
}

complex<double> XDecAmpBase::XToEtaFAmp(Spin lamX, Spin jf, EvtData* theData, 
					std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFMag, 
					std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFPhi,
					double fMass, double fWidth){

  complex<double> result(0.,0.);

  Vector4<double > p4PiPi(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
			  theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
			  theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
			  theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());  

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;
  
  for ( itXMag=XToEtaFMag.begin(); itXMag!=XToEtaFMag.end(); ++itXMag){
    boost::shared_ptr<const JPCLS> XState=itXMag->first;
    double theXMag=itXMag->second;
    double theXPhi=XToEtaFPhi[XState];
    complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
    complex<double> amp(0.,0.);     
    for(Spin lamf = -jf; lamf <= jf; lamf++){
      
      if( fabs(lamf)> _spinX || fabs(lamf)>XState->S) continue;
      
      amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)
	*Clebsch(XState->L, 0, XState->S, lamf, _spinX, lamf)
	*Clebsch(jf, lamf, 0, 0, XState->S, lamf)
	*conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_XTofEtadec][_spinX][lamX][lamf])
	*conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_fToPiPiDec][jf][lamf][0])
	*BreitWignerBlattW(p4PiPi, _massPi0, _massPi0, fMass, fWidth, jf);
     }
     
     result+= amp;
   }
   
   return result;

}


void XDecAmpBase::print(std::ostream& os) const{
  return; //dummy
}

void XDecAmpBase::initialize(){
  std::vector<std::string>::const_iterator it;

  for (it=_hypVec.begin(); it!=_hypVec.end(); ++it){

    if (it->compare(0, _a980piKey.size(), _a980piKey)== 0){
      Info << "hypothesis\t" << _a980piKey << "\t enabled" << endmsg;
      _a980piHyp=true;
    }

    else if (it->compare(0, _sigmaEtaKey.size(), _sigmaEtaKey)== 0){
      Info << "hypothesis\t" << _sigmaEtaKey << "\t enabled" << endmsg;
      _sigmaEtaHyp=true;
    }

  }

}
