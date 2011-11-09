#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsPsi2SToKpKmPiGamLh::AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList) :
  _Psi2SToKpKmPiGamEvtListPtr(theEvtList)
{
  _evtDataVec=_Psi2SToKpKmPiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2SToKpKmPiGamEvtListPtr->getMcVecs();
}

AbsPsi2SToKpKmPiGamLh::AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theAbsPsi2SToKpKmPiGamLhPtr):
  _Psi2SToKpKmPiGamEvtListPtr(theAbsPsi2SToKpKmPiGamLhPtr->getEventList())
{
  _evtDataVec=_Psi2SToKpKmPiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2SToKpKmPiGamEvtListPtr->getMcVecs();
}

AbsPsi2SToKpKmPiGamLh::~AbsPsi2SToKpKmPiGamLh()
{
}

double AbsPsi2SToKpKmPiGamLh::calcLogLh(fitParams& theParamVal){
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*>::iterator iterm;
  for (iterm=_evtMCVec.begin(); iterm!=_evtMCVec.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/_evtMCVec.size());

  logLH=_evtDataVec.size()/2.*(LH_mc/_evtMCVec.size()-1)*(LH_mc/_evtMCVec.size()-1)
    -logLH_data
    +_evtDataVec.size()*logLH_mc_Norm;

  Info << "current LH = " << logLH << endmsg;

 return logLH;

}

double AbsPsi2SToKpKmPiGamLh::calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, fitParams& theParamVal){

  double phaseSpaceVal=theParamVal.otherParams[paramEnumChic1ToKpKmPiGam::phaseSpace];  


  Spin Psi2SM=1;
  Spin GamM=1;
  complex<double> AmpPsi2SMpGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=1;
  GamM=-1;
  complex<double> AmpPsi2SMpGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=-1;
  GamM=1; 
  complex<double> AmpPsi2SMmGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=-1;
  GamM=-1; 
  complex<double> AmpPsi2SMmGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

//   DebugMsg << "AmpPsi2SMp " << AmpPsi2SMp << endmsg;
  
  double result=norm(AmpPsi2SMpGp)+norm(AmpPsi2SMpGm)+norm(AmpPsi2SMmGp)+norm(AmpPsi2SMmGm)+phaseSpaceVal;

  return result;  
}



complex<double> AbsPsi2SToKpKmPiGamLh::a980Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiPhi, double a980Mass, double a980Width, Spin& lamChi){

  Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());
  

  complex<double> result=conj(theData->DfChiToa0Pi[1][lamChi][0])*BreitWigner(KpKm4V, a980Mass, a980Width); // Amp(a0(980)->K+K-)=1;
  
  complex<double> ChiToA980pi0AmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToa0PiMag;
  
  for ( itChiToa0PiMag=ChiToa0PiMag.begin(); itChiToa0PiMag!=ChiToa0PiMag.end(); ++itChiToa0PiMag){
    boost::shared_ptr<const JPCLS> ChiToa0State=itChiToa0PiMag->first;
    double theChiToa0Mag=itChiToa0PiMag->second;
    double theChiToa0Phi=ChiToa0PiPhi[ChiToa0State];
    complex<double> expiphiChiToa0(cos(theChiToa0Phi), sin(theChiToa0Phi));
    
    ChiToA980pi0AmpTmpls+=theChiToa0Mag*expiphiChiToa0*sqrt(2*ChiToa0State->L+1); 
    //Clebsch(ChiToa0State->L=1,0,ChiToa0State->S=0, a0lam-pilam=0, ChiToa0State->J=1, 0)=1  
    //Clebsch(a0.J=0, a0lam=0, piJ=0, pilam=0, ChiToa0State->S=0, a0lam-pilam=0)=1
  }
  
  result*=ChiToA980pi0AmpTmpls;

  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::a980FlatteAmp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, fitParams& theParamVal, Spin& lamChi){

  Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToa0PiMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::a980Pi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToa0PiPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::a980Pi];  
  double a980FlatteMass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::a980];
  double a980FlatgKK=theParamVal.gFactors[paramEnumChic1ToKpKmPiGam::a980gKK];
  double a980FlatgEtaPi=theParamVal.gFactors[paramEnumChic1ToKpKmPiGam::a980gEtaPi];

  complex<double> result=conj(theData->DfChiToa0Pi[1][lamChi][0])*FlatteA980(KpKm4V, a980FlatteMass, a980FlatgKK, a980FlatgEtaPi); // Amp(a0(980)->K+K-)=1;  

  complex<double> ChiToA980pi0AmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToa0PiMag;
  
  for ( itChiToa0PiMag=ChiToa0PiMag.begin(); itChiToa0PiMag!=ChiToa0PiMag.end(); ++itChiToa0PiMag){
    boost::shared_ptr<const JPCLS> ChiToa0State=itChiToa0PiMag->first;
    double theChiToa0Mag=itChiToa0PiMag->second;
    double theChiToa0Phi=ChiToa0PiPhi[ChiToa0State];
    complex<double> expiphiChiToa0(cos(theChiToa0Phi), sin(theChiToa0Phi));
    
    ChiToA980pi0AmpTmpls+=theChiToa0Mag*expiphiChiToa0*sqrt(2*ChiToa0State->L+1); 
    //Clebsch(ChiToa0State->L=1,0,ChiToa0State->S=0, a0lam-pilam=0, ChiToa0State->J=1, 0)=1  
    //Clebsch(a0.J=0, a0lam=0, piJ=0, pilam=0, ChiToa0State->S=0, a0lam-pilam=0)=1
  }
  
  result*=ChiToA980pi0AmpTmpls;

  return result;
}


// complex<double> AbsPsi2SToKpKmPiGamLh::a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi){

//   complex<double> result(0.,0.);
//   Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());

//   for (Spin lam_a2=-2; lam_a2<=2; lam_a2++){
    
//     complex<double> ChiToa2PiAmpTmp1=conj(theData->DfChiToa2Pi[1][lamChi][lam_a2])
//       *sqrt(5)*conj(theData->Dfa2[2][lam_a2][0])*BreitWigner(KpKm4V, a2Mass, a2Width);
    
//     complex<double> ChiToa2PiAmpTmpls(0.,0.);
    
//     std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToa2Pi;
//     for ( itChiToa2Pi=ChiToa2Pi.begin(); itChiToa2Pi!=ChiToa2Pi.end(); ++itChiToa2Pi){
// 	    boost::shared_ptr<const JPCLS> ChiToa2PiState=itChiToa2Pi->first;
	    
// 	    if ( fabs(lam_a2)>ChiToa2PiState->J || fabs(lam_a2)>ChiToa2PiState->S) continue;
	    
// 	    double theChiToa2PiMag=itChiToa2Pi->second.first;
// 	    double theChiToa2PiPhi=itChiToa2Pi->second.second;
// 	    complex<double> expiphiChiToa2Pi(cos(theChiToa2PiPhi), sin(theChiToa2PiPhi));
	    
	    
// 	    ChiToa2PiAmpTmpls+=theChiToa2PiMag*expiphiChiToa2Pi*sqrt(2*ChiToa2PiState->L+1)
// 	      *Clebsch(ChiToa2PiState->L,0,ChiToa2PiState->S, lam_a2, ChiToa2PiState->J, lam_a2)
// 	      *Clebsch(2, lam_a2, 0, 0, ChiToa2PiState->S, lam_a2);
//           }
	  
// 	  ChiToa2PiAmpTmp1*=ChiToa2PiAmpTmpls;
// 	  result+=ChiToa2PiAmpTmp1;
//   }

//   return result;
// }



complex<double> AbsPsi2SToKpKmPiGamLh::K0_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KPhi , double K0_Mass, double K0_Width, Spin& lamChi){

  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());


  complex<double> result=conj(theData->DfChiToKpPi[1][lamChi][0])
    *BreitWigner(KpPi4V, K0_Mass, K0_Width)
    +conj(theData->DfChiToKmPi[1][lamChi][0])
    *BreitWigner(KmPi4V, K0_Mass, K0_Width);


  complex<double> ChiToK1400_0_KAmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToK1400_0_KMag;
  for ( itChiToK1400_0_KMag=ChiToK_0_KMag.begin(); itChiToK1400_0_KMag!=ChiToK_0_KMag.end(); ++itChiToK1400_0_KMag){
 
    boost::shared_ptr<const JPCLS> ChiToK1400_0_KState=itChiToK1400_0_KMag->first;
    double theChiToK1400_0_KMag=itChiToK1400_0_KMag->second;
    double theChiToK1400_0_KPhi=ChiToK_0_KPhi[ChiToK1400_0_KState];
    complex<double> expiphiChiToK1400_0_K(cos(theChiToK1400_0_KPhi), sin(theChiToK1400_0_KPhi));
    
    ChiToK1400_0_KAmpTmpls+=theChiToK1400_0_KMag*expiphiChiToK1400_0_K*sqrt(2*ChiToK1400_0_KState->L+1);
    //      *Clebsch(ChiToK1400_0_KState->L,0,ChiToK1400_0_KState->S, 0, ChiToK1400_0_KState->J, 0)
    //      *Clebsch(0, 0, 0, 0, 0, 0)
  }

  result*=ChiToK1400_0_KAmpTmpls;


  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::K1_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1_KPhi , double K1_Mass, double K1_Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

  for (Spin lamK1=-1; lamK1<=1; lamK1++){
    
//     complex<double> ChiToK1_1400KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1400])
//       *sqrt(3)*conj(theData->DfKp1[1][lamK1400][0])*BreitWigner(KpPi4V, K1_1400Mass, K1_1400Width)
//       +conj(theData->DfChiToKmPi[1][lamChi][lamK1400])
//       *sqrt(3)*conj(theData->DfKm1[1][lamK1400][0])*BreitWigner(KmPi4V, K1_1400Mass, K1_1400Width);

    complex<double> ChiToK1_KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1])
      *sqrt(3)*conj(theData->DfKp1[1][lamK1][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K1_Mass, K1_Width, 1)
      +conj(theData->DfChiToKmPi[1][lamChi][lamK1])
      *sqrt(3)*conj(theData->DfKm1[1][lamK1][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K1_Mass, K1_Width, 1);
    complex<double> ChiToK1_KAmpTmpls(0.,0.);

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToK1_KMag;
    for ( itChiToK1_KMag=ChiToK1_KMag.begin(); itChiToK1_KMag!=ChiToK1_KMag.end(); ++itChiToK1_KMag){
      boost::shared_ptr<const JPCLS> ChiToK1_KState=itChiToK1_KMag->first;
      double theChiToK1_KMag=itChiToK1_KMag->second;
      double theChiToK1_KPhi=ChiToK1_KPhi[ChiToK1_KState];
      complex<double> expiphiChiToK1_K(cos(theChiToK1_KPhi), sin(theChiToK1_KPhi));
      
      ChiToK1_KAmpTmpls+=theChiToK1_KMag*expiphiChiToK1_K*sqrt(2*ChiToK1_KState->L+1)
	*Clebsch(ChiToK1_KState->L,0,ChiToK1_KState->S, lamK1, ChiToK1_KState->J, lamK1)
	*Clebsch(1, lamK1, 0, 0, ChiToK1_KState->S, lamK1);
    }
    
    ChiToK1_KAmpTmp1*=ChiToK1_KAmpTmpls;
    result+=ChiToK1_KAmpTmp1;
  }

  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::K2_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KPhi , double K2_1400Mass, double K2_1400Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

  for (Spin lamK1400_2=-2; lamK1400_2<=2; lamK1400_2++){
    
    complex<double> ChiToK1400_2_KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1400_2])
      *sqrt(5)*conj(theData->DfKp2[2][lamK1400_2][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K2_1400Mass, K2_1400Width, 2)
      +conj(theData->DfChiToKmPi[1][lamChi][lamK1400_2])
      *sqrt(5)*conj(theData->DfKm2[2][lamK1400_2][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K2_1400Mass, K2_1400Width, 2);
    
    complex<double> ChiToK1400_2_KAmpTmpls(0.,0.);
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToK1400_2_KMag;
    for ( itChiToK1400_2_KMag=ChiToK1400_2_KMag.begin(); itChiToK1400_2_KMag!=ChiToK1400_2_KMag.end(); ++itChiToK1400_2_KMag){
	    boost::shared_ptr<const JPCLS> ChiToK1400_2_KState=itChiToK1400_2_KMag->first;
	    
	    if ( fabs(lamK1400_2)>ChiToK1400_2_KState->J || fabs(lamK1400_2)>ChiToK1400_2_KState->S) continue;
	    
	    double theChiToK1400_2_KMag=itChiToK1400_2_KMag->second;
	    double theChiToK1400_2_KPhi=ChiToK1400_2_KPhi[ChiToK1400_2_KState];
	    complex<double> expiphiChiToK1400_2_K(cos(theChiToK1400_2_KPhi), sin(theChiToK1400_2_KPhi));
	    
	    
	    ChiToK1400_2_KAmpTmpls+=theChiToK1400_2_KMag*expiphiChiToK1400_2_K*sqrt(2*ChiToK1400_2_KState->L+1)
	      *Clebsch(ChiToK1400_2_KState->L,0,ChiToK1400_2_KState->S, lamK1400_2, ChiToK1400_2_KState->J, lamK1400_2)
	      *Clebsch(2, lamK1400_2, 0, 0, ChiToK1400_2_KState->S, lamK1400_2);
          }
	  
	  ChiToK1400_2_KAmpTmp1*=ChiToK1400_2_KAmpTmpls;
	  result+=ChiToK1400_2_KAmpTmp1;
  }

  return result;
}


//   complex<double> AbsPsi2SToKpKmPiGamLh::KappaPole_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappa_K, double KappaMass, double KappaWidth, Spin& lamChi){

//     Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
//     Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

//     complex<double> i(0.,1.);    
    
//     complex<double> result=conj(theData->DfChiToKpPi[1][lamChi][0])
//       *KappaMass*KappaWidth/((KappaMass-i*KappaWidth)*(KappaMass-i*KappaWidth)-KpPi4V*KpPi4V)
//       +conj(theData->DfChiToKmPi[1][lamChi][0])
//       *KappaMass*KappaWidth/((KappaMass-i*KappaWidth)*(KappaMass-i*KappaWidth)-KmPi4V*KmPi4V);

//   complex<double> ChiToKappa_KAmpTmpls(0.,0.);
  
//   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToKappa_K;
//   for ( itChiToKappa_K=ChiToKappa_K.begin(); itChiToKappa_K!=ChiToKappa_K.end(); ++itChiToKappa_K){
//     boost::shared_ptr<const JPCLS> ChiToKappa_KState=itChiToKappa_K->first;
//     double theChiToKappa_KMag=itChiToKappa_K->second.first;
//     double theChiToKappa_KPhi=itChiToKappa_K->second.second;
//     complex<double> expiphiChiToKappa_K(cos(theChiToKappa_KPhi), sin(theChiToKappa_KPhi));
    
//     ChiToKappa_KAmpTmpls+=theChiToKappa_KMag*expiphiChiToKappa_K*sqrt(2*ChiToKappa_KState->L+1);
//     //      *Clebsch(ChiToK1400_0_KState->L,0,ChiToK1400_0_KState->S, 0, ChiToK1400_0_KState->J, 0)
//     //      *Clebsch(0, 0, 0, 0, 0, 0)
//   }

//   result*=ChiToKappa_KAmpTmpls;
//   return result;

// }



// complex<double> AbsPsi2SToKpKmPiGamLh::KKPi_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi, Spin& lamChi){

//   complex<double> result=conj(theData->DfChiToKKPi[1][lamChi][0]);   

//   complex<double> ChiToKKPiAmpTmpls(0.,0.);
  
//   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToKKPi;
  
//   for ( itChiToKKPi=ChiToKKPi.begin(); itChiToKKPi!=ChiToKKPi.end(); ++itChiToKKPi){
//     boost::shared_ptr<const JPCLS> ChiToKKPiState=itChiToKKPi->first;
//     double theChiToKKPiMag=itChiToKKPi->second.first;
//     double theChiToKKPiPhi=itChiToKKPi->second.second;
//     complex<double> expiphiChiToKKPi(cos(theChiToKKPiPhi), sin(theChiToKKPiPhi));
    
//     ChiToKKPiAmpTmpls+=theChiToKKPiMag*expiphiChiToKKPi*sqrt(2.*ChiToKKPiState->L+1); 
//     //Clebsch(ChiToKKPiState->L=1,0,ChiToKKPiState->S=0, lam=0, ChiToKKPiState->J=1, 0)=1  
//     //Clebsch(0, 0, 0, 0, ChiToKKPiState->S=0, 0)=1
//   }
  
//   result*=ChiToKKPiAmpTmpls;

//   return result;
// }

void AbsPsi2SToKpKmPiGamLh::print(std::ostream& os) const{
  os << "AbsPsi2SToKpKmPiGamLh::print\n";
}

