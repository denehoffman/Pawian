#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsPsi2STo2K2PiGamLh::AbsPsi2STo2K2PiGamLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList) :
  _Psi2STo2K2PiGamEvtListPtr(theEvtList),
  _fitParams2K2PiGam()
{
  _evtDataVec=_Psi2STo2K2PiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2STo2K2PiGamEvtListPtr->getMcVecs();
}

AbsPsi2STo2K2PiGamLh::AbsPsi2STo2K2PiGamLh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theAbsPsi2STo2K2PiGamLhPtr):
  _Psi2STo2K2PiGamEvtListPtr(theAbsPsi2STo2K2PiGamLhPtr->getEventList()),
  _fitParams2K2PiGam()
{
  _evtDataVec=_Psi2STo2K2PiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2STo2K2PiGamEvtListPtr->getMcVecs();
}

AbsPsi2STo2K2PiGamLh::~AbsPsi2STo2K2PiGamLh()
{
}

double AbsPsi2STo2K2PiGamLh::calcLogLh(const param2K2PiGam& theParamVal){
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*>::iterator iterm;
  for (iterm=_evtMCVec.begin(); iterm!=_evtMCVec.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log(LH_mc/_evtMCVec.size());

  logLH=_evtDataVec.size()*(LH_mc/_evtMCVec.size()-1)*(LH_mc/_evtMCVec.size()-1)
    -2.*logLH_data
    +2.*_evtDataVec.size()*logLH_mc_Norm;

  Info << "current LH = " << logLH << endmsg;

 return logLH;

}

double AbsPsi2STo2K2PiGamLh::calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const param2K2PiGam& theParamVal){
  double phaseSpaceVal=theParamVal.phaseSpace;
  
  complex<double> theDecAmp=chi0DecAmps(theParamVal, theData);

  Spin Psi2SM=1;
  Spin GamM=1;
  complex<double> AmpPsi2SMpGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

  Psi2SM=1;
  GamM=-1;
  complex<double> AmpPsi2SMpGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

  Psi2SM=-1;
  GamM=1; 
  complex<double> AmpPsi2SMmGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

  Psi2SM=-1;
  GamM=-1; 
  complex<double> AmpPsi2SMmGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

//   DebugMsg << "AmpPsi2SMp " << AmpPsi2SMp << endmsg;
  
  double result=norm(AmpPsi2SMpGp)+norm(AmpPsi2SMpGm)+norm(AmpPsi2SMmGp)+norm(AmpPsi2SMmGm)+phaseSpaceVal;

  return result;  
}
  
complex<double> AbsPsi2STo2K2PiGamLh::calcCoherentAmp(Spin Minit, Spin lamGam, const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;

  complex<double> PsiDecAmpTmpls(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
  for ( itPsi=PsiToChiGam.begin(); itPsi!=PsiToChiGam.end(); ++itPsi){
    
    boost::shared_ptr<const JPCLS> PsiState=itPsi->first;

    double thePsiMag=itPsi->second.first;
    double thePsiPhi=itPsi->second.second;
    complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));


    PsiDecAmpTmpls+=thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
      *Clebsch(PsiState->L,0,PsiState->S, -lamGam, PsiState->J, -lamGam)
      *Clebsch(0, 0, 1, -lamGam, PsiState->S, -lamGam);
  }

  complex<double> result=conj(theData->DfPsi[1][Minit][-lamGam])*PsiDecAmpTmpls;
  return result; 
}



complex<double> AbsPsi2STo2K2PiGamLh::chiTo2K892Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K892, double K892Mass, double K892Width){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  for (Spin lamK892_1=-1; lamK892_1<=1; lamK892_1++){

    Spin lamK892_2=lamK892_1;
    
    complex<double> ampKpPi0= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi0[1][lamK892_1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    complex<double> ampKpPi1= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi1[1][lamK892_1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    complex<double> ampKmPi0= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi0[1][lamK892_2][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    complex<double> ampKmPi1= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi1[1][lamK892_2][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiTo2K892;
    for ( itChiTo2K892=ChiTo2K892.begin(); itChiTo2K892!=ChiTo2K892.end(); ++itChiTo2K892){

      boost::shared_ptr<const JPCLS> ChiTo2K892State=itChiTo2K892->first;      
      double theMag=itChiTo2K892->second.first;
      double thePhi=itChiTo2K892->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));

      complex<double> ChiTo2K892Tmpls=ampKpPi0*ampKmPi1;
      ChiTo2K892Tmpls+=ampKpPi1*ampKmPi0;      

      ChiTo2K892Tmpls=theMag*expiphi*ChiTo2K892Tmpls*sqrt(2*ChiTo2K892State->L+1)
	 *Clebsch(ChiTo2K892State->L,0,ChiTo2K892State->S, 0, 0, 0)
	 *Clebsch(1, lamK892_1, 1, -lamK892_2, ChiTo2K892State->S, 0);

      result+=ChiTo2K892Tmpls;		
    }
  }

  return result;
}






complex<double> AbsPsi2STo2K2PiGamLh::chiToK1K1Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK1K1, double K1Mass0, double K1Width0,  double K1Mass1, double K1Width1){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  for (Spin lamK1=-1; lamK1<=1; lamK1++){

    complex<double> amp_K10_KpPi0= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi0[1][lamK1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KpPi0= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi0[1][lamK1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));    


    complex<double> amp_K10_KpPi1= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi1[1][lamK1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KpPi1= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi1[1][lamK1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));


    complex<double> amp_K10_KmPi0= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi0[1][lamK1][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KmPi0= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi0[1][lamK1][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));


    complex<double> amp_K10_KmPi1= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi1[1][lamK1][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KmPi1= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi1[1][lamK1][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    for ( it=ChiToK1K1.begin(); it!=ChiToK1K1.end(); ++it){

      boost::shared_ptr<const JPCLS> ChiToK1K1State=it->first;      
      double theMag=it->second.first;
      double thePhi=it->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));

      complex<double> ChiToK1K1Tmpls=amp_K10_KpPi0*amp_K11_KmPi1;
      ChiToK1K1Tmpls+=amp_K10_KpPi1*amp_K11_KmPi0;      
      ChiToK1K1Tmpls+=amp_K10_KmPi0*amp_K11_KpPi1;
      ChiToK1K1Tmpls+=amp_K10_KmPi1*amp_K11_KpPi0;

      result+=theMag*expiphi*ChiToK1K1Tmpls*sqrt(2*ChiToK1K1State->L+1)
	 *Clebsch(ChiToK1K1State->L,0,ChiToK1K1State->S, 0, 0, 0)
	 *Clebsch(1, lamK1, 1, -lamK1, ChiToK1K1State->S, 0);

    }
  }

  return result;
}









complex<double> AbsPsi2STo2K2PiGamLh::chiTo2K_2_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_2, double K_2_Mass, double K_2_Width){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  for (Spin lamK2_1=-2; lamK2_1<=2; lamK2_1++){

    Spin lamK2_2=lamK2_1;
    
    complex<double> ampKpPi0= 
      sqrt(5)*( conj(theData->DfKst2pToKpPi0[2][lamK2_1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    complex<double> ampKpPi1= 
      sqrt(5)*( conj(theData->DfKst2pToKpPi1[2][lamK2_1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    complex<double> ampKmPi0= 
      sqrt(5)*( conj(theData->DfKst2mToKmPi0[2][lamK2_2][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    complex<double> ampKmPi1= 
      sqrt(5)*( conj(theData->DfKst2mToKmPi1[2][lamK2_2][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiTo2K2;
    for ( itChiTo2K2=ChiTo2K_2.begin(); itChiTo2K2!=ChiTo2K_2.end(); ++itChiTo2K2){

      boost::shared_ptr<const JPCLS> ChiTo2K_2_State=itChiTo2K2->first;      
      double theMag=itChiTo2K2->second.first;
      double thePhi=itChiTo2K2->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));

      complex<double> ChiTo2K2Tmpls=ampKpPi0*ampKmPi1;
      ChiTo2K2Tmpls+=ampKpPi1*ampKmPi0;      

      ChiTo2K2Tmpls=theMag*expiphi*ChiTo2K2Tmpls*sqrt(2*ChiTo2K_2_State->L+1)
	*Clebsch(ChiTo2K_2_State->L,0,ChiTo2K_2_State->S, 0, 0, 0)
	*Clebsch(2, lamK2_1, 2, -lamK2_2, ChiTo2K_2_State->S, 0);

      result+=ChiTo2K2Tmpls;		
    }
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK0K2Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK0K2, double K0Mass, double K0Width,  double K2Mass, double K2Width){
  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

    complex<double> ampK0ToKpPi0=  sqrt(3)*BreitWigner(KpPi0, K0Mass, K0Width);    
    complex<double> ampK0ToKpPi1=  sqrt(3)*BreitWigner(KpPi1, K0Mass, K0Width);
    complex<double> ampK0ToKmPi0=  sqrt(3)*BreitWigner(KmPi0, K0Mass, K0Width);    
    complex<double> ampK0ToKmPi1=  sqrt(3)*BreitWigner(KmPi1, K0Mass, K0Width);
    complex<double> ampK2ToKpPi0=  sqrt(5)*(conj(theData->DfKst2pToKpPi0[2][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K2Mass, K2Width, 2));
    complex<double> ampK2ToKpPi1=  sqrt(5)*(conj(theData->DfKst2pToKpPi1[2][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K2Mass, K2Width, 2));
    complex<double> ampK2ToKmPi0=  sqrt(5)*(conj(theData->DfKst2mToKmPi0[2][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K2Mass, K2Width, 2));
    complex<double> ampK2ToKmPi1=  sqrt(5)*(conj(theData->DfKst2mToKmPi1[2][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K2Mass, K2Width, 2));

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToK0K2.begin(); it!=ChiToK0K2.end(); ++it){
    
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));



    //    result+=theMag*expiphi*
    //  (ampK0ToKpPi0*ampK2ToKmPi1+ampK0ToKpPi1*ampK2ToKmPi0+ampK0ToKmPi0+ampK2ToKpPi1+ampK0ToKmPi1*ampK2ToKpPi0);
    
    result+=theMag*expiphi*
      (ampK0ToKpPi0*ampK2ToKmPi1+ampK0ToKpPi1*ampK2ToKmPi0+ampK0ToKmPi0*ampK2ToKpPi1+ampK0ToKmPi1*ampK2ToKpPi0); 

  } 
 
  return result;

}

complex<double> AbsPsi2STo2K2PiGamLh::chiTo2K_0_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_0, double K_0_Mass0, double K_0_Width0, double K_0_Mass1, double K_0_Width1){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiTo2K0;
  for ( itChiTo2K0=ChiTo2K_0.begin(); itChiTo2K0!=ChiTo2K_0.end(); ++itChiTo2K0){
    
    double theMag=itChiTo2K0->second.first;
    double thePhi=itChiTo2K0->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));
    
    result+=theMag*expiphi*
      ( BreitWigner(KpPi0, K_0_Mass0, K_0_Width0) * BreitWigner(KmPi1, K_0_Mass1, K_0_Width1)
	+ BreitWigner(KpPi0, K_0_Mass1, K_0_Width1) * BreitWigner(KmPi1, K_0_Mass0, K_0_Width0)
	+ BreitWigner(KpPi1, K_0_Mass0, K_0_Width0) * BreitWigner(KmPi0, K_0_Mass1, K_0_Width1)
	+ BreitWigner(KpPi1, K_0_Mass1, K_0_Width1) * BreitWigner(KmPi0, K_0_Mass0, K_0_Width0) );
    
  } 
 
  return result;
}


// complex<double> AbsPsi2STo2K2PiGamLh::chiToK1ToK1piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1ToK1dPi, double K1Mass, double K1Width, double K1dMass, double K1dWidth){
//   complex<double> result(0.,0.);
  
//   Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
// 			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

//   Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
// 			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

//   Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
// 			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

//   Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
// 			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

//   Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
// 			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

//   Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
// 			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());

//   Spin lamK1=0;
//   Spin lamK1d=0;


//   complex<double> resultTmpLoop(0.,0.);
  
//   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itK1ToK1dPi;
//   for ( itK1ToK1dPi=K1ToK1dPi.begin(); itK1ToK1dPi!=K1ToK1dPi.end(); ++itK1ToK1dPi){
    
//     boost::shared_ptr<const JPCLS> K1ToK1PiState=itK1ToK1dPi->first;
//     double theMag=itK1ToK1dPi->second.first;
//     double thePhi=itK1ToK1dPi->second.second;
//     complex<double> expiphi(cos(thePhi), sin(thePhi));

//     double K1prodClebsch=sqrt(5.);         
//     double K1decClebsch=sqrt(2.*K1ToK1PiState->L+1)*Clebsch(K1ToK1PiState->L,0, K1ToK1PiState->S, 0, K1ToK1PiState->J, 0)
//       *sqrt(5.);  //K1d decay

//     complex<double> ampKp1ToKp1dPi0=
//       conj(theData->DfK1400pToKstpPi0[1][0][0])*BreitWignerBlattW(KpPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L)//K1dec
//       *conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec

//     complex<double> ampKp1ToKp1dPi1=
//       conj(theData->DfK1400pToKstpPi1[1][0][0])*BreitWignerBlattW(KpPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L) //K1dec
//       *conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
    
//     complex<double> ampKm1ToKm1dPi0=
//       conj(theData->DfK1400mToKstmPi0[1][0][0])*BreitWignerBlattW(KmPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L) //K1dec
//       *conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec
    
//     complex<double> ampKm1ToKm1dPi1=
//       conj(theData->DfK1400mToKstmPi1[1][0][0])*BreitWignerBlattW(KmPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L) //K1dec
//       *conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

// //     complex<double> ampKp1ToKp1dPi0=
// //       conj(theData->DfK1400pToKstpPi0[1][0][0])*BreitWigner(KpPiPi, K1Mass, K1Width)//K1dec
// //       *conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec

// //     complex<double> ampKp1ToKp1dPi1=
// //       conj(theData->DfK1400pToKstpPi1[1][0][0])*BreitWigner(KpPiPi, K1Mass, K1Width) //K1dec
// //       *conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
    
// //     complex<double> ampKm1ToKm1dPi0=
// //       conj(theData->DfK1400mToKstmPi0[1][0][0])*BreitWigner(KmPiPi, K1Mass, K1Width) //K1dec
// //       *conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec
    
// //     complex<double> ampKm1ToKm1dPi1=
// //       conj(theData->DfK1400mToKstmPi1[1][0][0])*BreitWigner(KmPiPi, K1Mass, K1Width) //K1dec
// //       *conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

//     resultTmpLoop+=theMag*expiphi*K1prodClebsch*K1decClebsch*(ampKp1ToKp1dPi0+ampKp1ToKp1dPi1+ampKm1ToKm1dPi0+ampKm1ToKm1dPi1);
//   }
//   result=resultTmpLoop;
//   return result;
// }

complex<double> AbsPsi2STo2K2PiGamLh::chiToK1ToK1piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Dec, double K1Mass, double K1Width, double K1dMass, double K1dWidth){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());

  Spin lamK1=0;
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;

  for ( it=K1Prod.begin(); it!=K1Prod.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMagProd=it->second.first;
    double thePhiProd=it->second.second;
    complex<double> expiphiProd(cos(thePhiProd), sin(thePhiProd));

    complex<double> currentResultDec(0.,0.);

    for ( itDec=K1Dec.begin(); itDec!=K1Dec.end(); ++itDec){

      boost::shared_ptr<const JPCLS> theDecState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));

      complex<double> ampKp1ToKp1dPi0(0.,0.);
      complex<double> ampKp1ToKp1dPi1(0.,0.);
      complex<double> ampKm1ToKm1dPi0(0.,0.); 
      complex<double> ampKm1ToKm1dPi1(0.,0.);

      for (Spin lamKdec=-1; lamKdec<=1; ++lamKdec){
	double tmpVal=sqrt(3.)*Clebsch(theDecState->L, 0., theDecState->S, lamKdec, theDecState->J, lamKdec)*Clebsch(1,lamKdec, 0, 0, theDecState->S, lamKdec);	
	ampKp1ToKp1dPi0+=tmpVal*conj(theData->DfK1400pToKstpPi0[1][0][lamKdec])*conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][lamKdec][0]);
	ampKp1ToKp1dPi1+=tmpVal*conj(theData->DfK1400pToKstpPi1[1][0][lamKdec])*conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][lamKdec][0]);
	ampKm1ToKm1dPi0+=tmpVal*conj(theData->DfK1400mToKstmPi0[1][0][lamKdec])*conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][lamKdec][0]);
	ampKm1ToKm1dPi1+=tmpVal*conj(theData->DfK1400mToKstmPi1[1][0][lamKdec])*conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][lamKdec][0]);	
      }

      ampKp1ToKp1dPi0*=BreitWignerBlattW(KpPiPi, K1dMass, 0.1349766, K1Mass, K1Width, theDecState->L)*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
      ampKp1ToKp1dPi1*= BreitWignerBlattW(KpPiPi, K1dMass, 0.1349766, K1Mass, K1Width, theDecState->L)*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
      ampKm1ToKm1dPi0*=BreitWignerBlattW(KmPiPi, K1dMass, 0.1349766, K1Mass, K1Width, theDecState->L)*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
      ampKm1ToKm1dPi1*=BreitWignerBlattW(KmPiPi, K1dMass, 0.1349766, K1Mass, K1Width, theDecState->L)*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

      currentResultDec+=theMagDec*expiphiDec*sqrt(2.*theState->L+1.)*Clebsch(theState->L, 0., theState->S, 0, theState->J, 0)*Clebsch(1,0, 0, 0, theState->S, 0)*(ampKp1ToKp1dPi0+ampKp1ToKp1dPi1+ampKm1ToKm1dPi0+ampKm1ToKm1dPi1);
    }

    result+=theMagProd*expiphiProd*currentResultDec;
  }

  return result;
}


complex<double> AbsPsi2STo2K2PiGamLh::chiToK1ToK0piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Dec, double K1Mass, double K1Width, double K0Mass, double K0Width){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itProd;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;
  for ( itProd=K1Prod.begin(); itProd!=K1Prod.end(); ++itProd){
    
    boost::shared_ptr<const JPCLS> prodState=itProd->first;
    double theMag=itProd->second.first;
    double thePhi=itProd->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K1prodClebsch=sqrt(3.);         
    double K1decClebsch=sqrt(3.);  //K1d decay   


    for ( itDec=K1Dec.begin(); itDec!=K1Dec.end(); ++itDec){
    
      boost::shared_ptr<const JPCLS> decState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));
      
       
      complex<double> ampKpToKp0Pi0=
	conj(theData->DfK1400pToKstpPi0[1][0][0])*BreitWignerBlattW(KpPiPi, K0Mass, 0.1349766, K1Mass, K1Width, decState->L)//K1dec
	*BreitWigner(KpPi1, K0Mass, K0Width); 
      
      complex<double> ampKp1400ToKp1dPi1=
	conj(theData->DfK1400pToKstpPi1[1][0][0])*BreitWignerBlattW(KpPiPi, K0Mass, 0.1349766, K1Mass, K1Width, decState->L) //K1dec
	*BreitWigner(KpPi0, K0Mass, K0Width);
      
      complex<double> ampKm1400ToKm1dPi0=
	conj(theData->DfK1400mToKstmPi0[1][0][0])*BreitWignerBlattW(KmPiPi, K0Mass, 0.1349766, K1Mass, K1Width, decState->L) //K1dec
	*BreitWigner(KmPi1, K0Mass, K0Width); //K1ddec
      
      complex<double> ampKm1400ToKp1dPi1=
	conj(theData->DfK1400mToKstmPi1[1][0][0])*BreitWignerBlattW(KmPiPi, K0Mass, 0.1349766, K1Mass, K1Width, decState->L) //K1dec
	*BreitWigner(KmPi0, K0Mass, K0Width); //K1ddec
      
      result+=theMag*expiphi*theMagDec*expiphiDec*K1prodClebsch*K1decClebsch*(ampKpToKp0Pi0+ampKp1400ToKp1dPi1+ampKm1400ToKm1dPi0+ampKm1400ToKp1dPi1);
    }
  }
  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK0KT0KpipiKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K0dMass, double K0dWidth){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=K0ToK0Pi.begin(); it!=K0ToK0Pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> K0ToK0PiState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K1prodClebsch=sqrt(3.);         
    double K1decClebsch=sqrt(3.);  //K1d decay   

    complex<double> ampKpToKp0Pi0=BreitWigner(KpPiPi, K0pMass, K0pWidth)//K0prod
      *BreitWigner(KpPi1, K0dMass, K0dWidth); 

    complex<double> ampKpToKp1dPi1=BreitWigner(KpPiPi, K0pMass, K0pWidth) //K0prod
      *BreitWigner(KpPi0, K0dMass, K0dWidth);

    complex<double> ampKmToKm0Pi0=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *BreitWigner(KmPi1, K0dMass, K0dWidth);     

    complex<double> ampKmToKm1Pi1=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *BreitWigner(KmPi0, K0dMass, K0dWidth);

    result+=theMag*expiphi*K1prodClebsch*K1decClebsch*(ampKpToKp0Pi0+ampKpToKp1dPi1+ampKmToKm0Pi0+ampKmToKm1Pi1);
  }

  return result;
}


complex<double> AbsPsi2STo2K2PiGamLh::chiToK0KpiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK0KPi, double K0Mass, double K0Width){

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToK0KPi.begin(); it!=ChiToK0KPi.end(); ++it){
    
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    complex<double> ampKpToKpPi0=BreitWigner(KpPi0, K0Mass, K0Width); 
    complex<double> ampKpToKpPi1=BreitWigner(KpPi1, K0Mass, K0Width); 
    complex<double> ampKmToKmPi0=BreitWigner(KmPi0, K0Mass, K0Width); 
    complex<double> ampKmToKmPi1=BreitWigner(KmPi1, K0Mass, K0Width); 


    result+=theMag*expiphi*(ampKpToKpPi0+ampKpToKpPi1+ampKmToKmPi0+ampKmToKmPi1);
  }

  return result;
}





complex<double> AbsPsi2STo2K2PiGamLh::chiToKst1KpiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToKst1KPi, double Kst1Mass, double Kst1Width){

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToKst1KPi.begin(); it!=ChiToKst1KPi.end(); ++it){
    
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    complex<double> ampKpToKpPi0=BreitWigner(KpPi0, Kst1Mass, Kst1Width)*conj(theData->DfKst1pToKpPi0[1][0][0]); 
    complex<double> ampKpToKpPi1=BreitWigner(KpPi1, Kst1Mass, Kst1Width)*conj(theData->DfKst1pToKpPi1[1][0][0]); 
    complex<double> ampKmToKmPi0=BreitWigner(KmPi0, Kst1Mass, Kst1Width)*conj(theData->DfKst1mToKmPi0[1][0][0]); 
    complex<double> ampKmToKmPi1=BreitWigner(KmPi1, Kst1Mass, Kst1Width)*conj(theData->DfKst1mToKmPi1[1][0][0]); 


    result+=theMag*expiphi*3.*(ampKpToKpPi0+ampKpToKpPi1+ampKmToKmPi0+ampKmToKmPi1);
  }

  return result;
}






complex<double> AbsPsi2STo2K2PiGamLh::chiToK0KT0K1piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K1dMass, double K1dWidth){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=K0ToK0Pi.begin(); it!=K0ToK0Pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> K0ToK0PiState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K0proddecClebsch=2.*sqrt(3.);         
    double K1decClebsch=sqrt(5.);  //K1d decay   

    complex<double> ampKpToK1pPi0=BreitWigner(KpPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); 


    complex<double> ampKpToK1pPi1=BreitWigner(KpPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);


    complex<double> ampKmToK1mPi0=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

    complex<double> ampKmToK1mPi1=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);


    result+=theMag*expiphi*K0proddecClebsch*K1decClebsch*(ampKpToK1pPi0+ampKpToK1pPi1+ampKmToK1mPi0+ampKmToK1mPi1);
  }

  return result;
}



complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif0_k, double f980_piMass, double f980_piWidth, double f0_kMass, double f0_kWidth){

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
			 theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());

  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
			 theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_pif0_k.begin(); it!=ChiTof980_pif0_k.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(PiPi, f980_piMass,f980_piWidth)*BreitWigner(KK, f0_kMass, f0_kWidth);
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif0_k, double f980_piMass, double f980_g1, double f980_g2, double f0_kMass, double f0_kWidth){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_pif0_k.begin(); it!=ChiTof980_pif0_k.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_piMass, f980_g1, f980_g2)*BreitWigner(KK, f0_kMass, f0_kWidth);
  }

  return result;

}

complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_kf0_piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_kf0_pi, double f980_Mass, double f980_gKK, double f980_gPiPi, double f0_piMass, double f0_piWidth){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_kf0_pi.begin(); it!=ChiTof980_kf0_pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(KK, decPair2, decPair1, f980_Mass, f980_gKK, f980_gPiPi)*BreitWigner(PiPi, f0_piMass, f0_piWidth);
  }

  return result;

}

complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_pif2_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif2_k, double f980_piMass, double f980_g1, double f980_g2, double f2_kMass, double f2_kWidth){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_pif2_k.begin(); it!=ChiTof980_pif2_k.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_piMass, f980_g1, f980_g2)*BreitWignerBlattW(KK, 0.493677, 0.493677, f2_kMass, f2_kWidth, 2)*conj(theData->Dff2ToKK[2][0][0]);
  }

  return result;


}

complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_kf2_piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_kf2_pi, double f980_Mass, double f980_gKK, double f980_gPiPi, double f2_piMass, double f2_piWidth){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_kf2_pi.begin(); it!=ChiTof980_kf2_pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(KK, decPair2, decPair1, f980_Mass, f980_gKK, f980_gPiPi)*BreitWignerBlattW(PiPi, 0.1349766, 0.1349766, f2_piMass, f2_piWidth, 2)*conj(theData->Dff2ToPiPi[2][0][0]);
  }

  return result;

}


complex<double> AbsPsi2STo2K2PiGamLh::chiTof980f980Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980980, double f980_piMass, double f980_g1, double f980_g2){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980980.begin(); it!=ChiTof980980.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_piMass, f980_g1, f980_g2)*Flatte(KK, decPair2, decPair1, f980_piMass, f980_g2, f980_g1);
  }

  return result;

}

complex<double>  AbsPsi2STo2K2PiGamLh::chiTof0_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof0_pif0_k, double f0_pi_Mass, double f0_pi_Width,  double f0_kMass, double f0_kWidth){
 
  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof0_pif0_k.begin(); it!=ChiTof0_pif0_k.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(PiPi, f0_pi_Mass, f0_pi_Width)*BreitWigner(KK, f0_kMass, f0_kWidth);
  }

  return result;
}

complex<double>  AbsPsi2STo2K2PiGamLh::chiTof0_pif2_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof0_pif2_k, double f0_pi_Mass, double f0_pi_Width, double f2_kMass, double f2_kWidth){
 
  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof0_pif2_k.begin(); it!=ChiTof0_pif2_k.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(PiPi, f0_pi_Mass, f0_pi_Width)*BreitWignerBlattW(KK, 0.493677, 0.493677, f2_kMass, f2_kWidth, 2)*conj(theData->Dff2ToKK[2][0][0]);
  }

  return result;
}

complex<double>  AbsPsi2STo2K2PiGamLh::chiTof2_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof2_pif0_k, double f2_pi_Mass, double f2_pi_Width, double f0_kMass, double f0_kWidth){
 
  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof2_pif0_k.begin(); it!=ChiTof2_pif0_k.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(KK, f0_kMass, f0_kWidth)*BreitWignerBlattW(PiPi, 0.1349766, 0.1349766, f2_pi_Mass, f2_pi_Width, 2)*conj(theData->Dff2ToPiPi[2][0][0]);
  }

  return result;
}


complex<double> AbsPsi2STo2K2PiGamLh::chiToK_0_KToKf980KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_0_KToKf980K, double f980_Mass, double f980_gPiPi,  double f980_gKK, double K_0_Mass, double K_0_Width){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);
  
  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  
  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);
  
  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());


  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToK_0_KToKf980K.begin(); it!=ChiToK_0_KToKf980K.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_Mass, f980_gPiPi, f980_gKK)*(BreitWigner(KpPiPi, K_0_Mass, K_0_Width)+BreitWigner(KmPiPi, K_0_Mass, K_0_Width));
  }

  return result;

}


complex<double> AbsPsi2STo2K2PiGamLh::chiToK_0_KToKf0KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_0_KToKf0K, double f0_Mass, double f0_Width, double K_0_Mass, double K_0_Width){

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());


  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToK_0_KToKf0K.begin(); it!=ChiToK_0_KToKf0K.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(PiPi, f0_Mass, f0_Width)*(BreitWigner(KpPiPi, K_0_Mass, K_0_Width)+BreitWigner(KmPiPi, K_0_Mass, K_0_Width));
  }

  return result;

}


complex<double> AbsPsi2STo2K2PiGamLh::chiToPi2Pi0Tof2PiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& Pi_2_Dec, double Pi_2_Mass, double Pi_2_Width, double f_2_Mass, double f_2_Width){

  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
		     theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());

  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
		     theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());

  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;

  for ( it=ChiToPi_2_Prod.begin(); it!=ChiToPi_2_Prod.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMagProd=it->second.first;
    double thePhiProd=it->second.second;
    complex<double> expiphiProd(cos(thePhiProd), sin(thePhiProd));

    complex<double> currentResultDec(0.,0.);
    for ( itDec=Pi_2_Dec.begin(); itDec!=Pi_2_Dec.end(); ++itDec){
      Spin lamPi2=0;      
      boost::shared_ptr<const JPCLS> theDecState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));
      complex<double> tmpDec(0.,0.);
      for (Spin lamf2=-2; lamf2<=2; ++lamf2){
	tmpDec+=sqrt(2.*theDecState->L+1)*Clebsch(theDecState->L, 0., theDecState->S, lamf2, theDecState->J, lamf2)*Clebsch(2,lamf2, 0, 0, theDecState->S, lamf2)
	  *sqrt(5.)
	  *(BreitWigner(KKPi0, Pi_2_Mass, Pi_2_Width)*theData->DfPi2Tof2Pi0[theDecState->J][lamPi2][lamf2]*theData->Dff2ToKKviaKKpi0[2][lamf2][0]
	    + BreitWigner(KKPi1, Pi_2_Mass, Pi_2_Width)*theData->DfPi2Tof2Pi1[theDecState->J][lamPi2][lamf2]*theData->Dff2ToKKviaKKpi1[2][lamf2][0]);      
      }
      currentResultDec+=tmpDec*theMagDec*expiphiDec;
    }

    result+=theMagProd*expiphiProd*currentResultDec*sqrt(2.*theState->L+1)*Clebsch(theState->L, 0., theState->S, 0, theState->J, 0);
  }

  result*=BreitWignerBlattW(KK, 0.493677, 0.493677, f_2_Mass, f_2_Width, 2);
  return result;
}

complex<double>  AbsPsi2STo2K2PiGamLh::chiToPi0PiToKstarKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_0_Dec, double Pi_0_Mass, double Pi_0_Width, double Kstar_Mass, double Kstar_Width){

  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
		     theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());

  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
		     theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());

  Vector4<double> Kstarp_pi0(theData->KpPi0_HeliKKPi0_4V.E(), theData->KpPi0_HeliKKPi0_4V.Px(), 
		     theData->KpPi0_HeliKKPi0_4V.Py(), theData->KpPi0_HeliKKPi0_4V.Pz());

  Vector4<double> Kstarp_pi1(theData->KpPi1_HeliKKPi1_4V.E(), theData->KpPi1_HeliKKPi1_4V.Px(), 
		     theData->KpPi1_HeliKKPi1_4V.Py(), theData->KpPi1_HeliKKPi1_4V.Pz());

  Vector4<double> Kstarm_pi0(theData->KmPi0_HeliKKPi0_4V.E(), theData->KmPi0_HeliKKPi0_4V.Px(), 
		     theData->KmPi0_HeliKKPi0_4V.Py(), theData->KmPi0_HeliKKPi0_4V.Pz());

  Vector4<double> Kstarm_pi1(theData->KmPi1_HeliKKPi1_4V.E(), theData->KmPi1_HeliKKPi1_4V.Px(), 
		     theData->KmPi1_HeliKKPi1_4V.Py(), theData->KmPi1_HeliKKPi1_4V.Pz());

  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;

  for ( itDec=ChiToPi_0_Dec.begin(); itDec!=ChiToPi_0_Dec.end(); ++itDec){

    boost::shared_ptr<const JPCLS> theState=itDec->first;
    double theMagDec=itDec->second.first;
    double thePhiDec=itDec->second.second;
    complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));

    complex<double> currentResultDec(0.,0.);

    complex<double> tmpDec(0.,0.);
    
    tmpDec+=sqrt(2.*theState->L+1)*Clebsch(theState->L, 0., theState->S, 0, theState->J, 0)*Clebsch(1,0, 0, 0, theState->S, 0)
      *(( BreitWigner(KKPi0, Pi_0_Mass, Pi_0_Width)
	  *( BreitWignerBlattW(Kstarp_pi0, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,1)
	     *theData->DfKst1pToKpPi0ViaKKPi0[1][0][0]
	     +BreitWignerBlattW(Kstarm_pi0, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,1)
	     *theData->DfKst1mToKmPi0ViaKKPi0[1][0][0]))
	+( BreitWigner(KKPi1, Pi_0_Mass, Pi_0_Width)
	   *( BreitWignerBlattW(Kstarp_pi1, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,1)
	      *theData->DfKst1pToKpPi1ViaKKPi1[1][0][0] 
	      +BreitWignerBlattW(Kstarm_pi1, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,1)
	      *theData->DfKst1mToKmPi1ViaKKPi1[1][0][0] ))
	);
    
    
    currentResultDec+=tmpDec*theMagDec*expiphiDec;
    
    result+=theMagDec*expiphiDec*currentResultDec;
  }

  return result;
}


complex<double>  AbsPsi2STo2K2PiGamLh::chiToPi2Pi0ToKstarKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& Pi_2_Dec, double Pi_2_Mass, double Pi_2_Width, double Kstar_Mass, double Kstar_Width, Spin Kstar_Spin){

  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
		     theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());

  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
		     theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());

  Vector4<double> Kstarp_pi0(theData->KpPi0_HeliKKPi0_4V.E(), theData->KpPi0_HeliKKPi0_4V.Px(), 
		     theData->KpPi0_HeliKKPi0_4V.Py(), theData->KpPi0_HeliKKPi0_4V.Pz());

  Vector4<double> Kstarp_pi1(theData->KpPi1_HeliKKPi1_4V.E(), theData->KpPi1_HeliKKPi1_4V.Px(), 
		     theData->KpPi1_HeliKKPi1_4V.Py(), theData->KpPi1_HeliKKPi1_4V.Pz());

  Vector4<double> Kstarm_pi0(theData->KmPi0_HeliKKPi0_4V.E(), theData->KmPi0_HeliKKPi0_4V.Px(), 
		     theData->KmPi0_HeliKKPi0_4V.Py(), theData->KmPi0_HeliKKPi0_4V.Pz());

  Vector4<double> Kstarm_pi1(theData->KmPi1_HeliKKPi1_4V.E(), theData->KmPi1_HeliKKPi1_4V.Px(), 
		     theData->KmPi1_HeliKKPi1_4V.Py(), theData->KmPi1_HeliKKPi1_4V.Pz());

  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;
  
  int L_KstarDec=int(Kstar_Spin);

  for ( it=ChiToPi_2_Prod.begin(); it!=ChiToPi_2_Prod.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMagProd=it->second.first;
    double thePhiProd=it->second.second;
    complex<double> expiphiProd(cos(thePhiProd), sin(thePhiProd));

    complex<double> currentResultDec(0.,0.);
    for ( itDec=Pi_2_Dec.begin(); itDec!=Pi_2_Dec.end(); ++itDec){
      Spin lamPi2=0;      
      boost::shared_ptr<const JPCLS> theDecState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));
      complex<double> tmpDec(0.,0.);

      for (Spin lamKstar=-Kstar_Spin; lamKstar<=Kstar_Spin; ++lamKstar){
      
	if( fabs(lamKstar)>theDecState->J || fabs(lamKstar)>theDecState->S) continue;
	
	tmpDec+=sqrt(2.*theDecState->L+1)*Clebsch(theDecState->L, 0, theDecState->S, lamKstar, theDecState->J, lamKstar)*Clebsch(Kstar_Spin,lamKstar, 0, 0, theDecState->S, lamKstar)*sqrt(2.*L_KstarDec+1.)
	  *(( BreitWigner(KKPi0, Pi_2_Mass, Pi_2_Width)
	      *( BreitWignerBlattW(Kstarp_pi0, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,L_KstarDec)
		 *theData->DfPi2ToKstarpK_pi0[theDecState->J][lamPi2][lamKstar]*theData->DfKst1pToKpPi0ViaKKPi0[Kstar_Spin][lamKstar][0]
		 +BreitWignerBlattW(Kstarm_pi0, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,L_KstarDec)
		 *theData->DfPi2ToKstarmK_pi0[theDecState->J][lamPi2][lamKstar]*theData->DfKst1mToKmPi0ViaKKPi0[Kstar_Spin][lamKstar][0]))
	    +( BreitWigner(KKPi1, Pi_2_Mass, Pi_2_Width)
	       *( BreitWignerBlattW(Kstarp_pi1, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,L_KstarDec)
		  *theData->DfPi2ToKstarpK_pi1[theDecState->J][lamPi2][lamKstar]*theData->DfKst1pToKpPi1ViaKKPi1[Kstar_Spin][lamKstar][0] 
		  +BreitWignerBlattW(Kstarm_pi1, 0.493677, 0.1349766, Kstar_Mass, Kstar_Width,L_KstarDec)
		  *theData->DfPi2ToKstarmK_pi1[theDecState->J][lamPi2][lamKstar]*theData->DfKst1mToKmPi1ViaKKPi1[Kstar_Spin][lamKstar][0] ))
	    );

      }
      currentResultDec+=tmpDec*theMagDec*expiphiDec*sqrt(2.*theState->L+1)*Clebsch(theState->L, 0., theState->S, 0, theState->J, 0);
    }

    result+=theMagProd*expiphiProd*currentResultDec;
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToPi0Pi0Tof980kkAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi0Pi0Amp, double f980_Mass, double f980_gKK, double f980_gPiPi, double Pi_0_Mass, double Pi_0_Width){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);
  
  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  
  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);
  
  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
		     theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());

  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
		     theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());

  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToPi0Pi0Amp.begin(); it!=ChiToPi0Pi0Amp.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(KK, decPair2, decPair1, f980_Mass, f980_gPiPi, f980_gKK)
      *(BreitWigner(KKPi0, Pi_0_Mass, Pi_0_Width)+BreitWigner(KKPi1, Pi_0_Mass, Pi_0_Width));
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToPi0Pi0Tof0Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi0Pi0Amp, double Pi_0_Mass, double Pi_0_Width, double f0Mass, double f0Width){

  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
		     theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());

  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
		     theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());

  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToPi0Pi0Amp.begin(); it!=ChiToPi0Pi0Amp.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));
    
    result+=theMag*expiphi*BreitWigner(KK, f0Mass, f0Width)
      *(BreitWigner(KKPi0, Pi_0_Mass, Pi_0_Width)+BreitWigner(KKPi1, Pi_0_Mass, Pi_0_Width));
  }

  return result;
}


complex<double> AbsPsi2STo2K2PiGamLh::chiToPi0Pi0ToK0KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi0Pi0Amp, double Pi_0_Mass, double Pi_0_Width, double K0Mass, double K0Width){

  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
		     theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());

  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
		     theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
		     theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
		     theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
		     theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
		     theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToPi0Pi0Amp.begin(); it!=ChiToPi0Pi0Amp.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));
    
    result+=theMag*expiphi
      *( BreitWigner(KKPi0, Pi_0_Mass, Pi_0_Width)
	 *(BreitWigner(KpPi0, K0Mass, K0Width)+BreitWigner(KmPi0, K0Mass, K0Width) ) 
	 +
	 BreitWigner(KKPi1, Pi_0_Mass, Pi_0_Width)
	 *(BreitWigner(KpPi1, K0Mass, K0Width)+BreitWigner(KmPi1, K0Mass, K0Width) )
	 );
  }

  return result;
}



complex<double> AbsPsi2STo2K2PiGamLh::chiToPi2Pi0Tof0PiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToPi_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& Pi_2_Dec, double Pi_2_Mass, double Pi_2_Width, double f0Mass, double f0Width){

  Vector4<double> KKPi0(theData->KKPi0_HeliChic0_4V.E(), theData->KKPi0_HeliChic0_4V.Px(), 
			theData->KKPi0_HeliChic0_4V.Py(), theData->KKPi0_HeliChic0_4V.Pz());
  
  Vector4<double> KKPi1(theData->KKPi1_HeliChic0_4V.E(), theData->KKPi1_HeliChic0_4V.Px(), 
			theData->KKPi1_HeliChic0_4V.Py(), theData->KKPi1_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;

  for ( it=ChiToPi_2_Prod.begin(); it!=ChiToPi_2_Prod.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMagProd=it->second.first;
    double thePhiProd=it->second.second;
    complex<double> expiphiProd(cos(thePhiProd), sin(thePhiProd));

    complex<double> currentResultDec(0.,0.);

    for ( itDec=Pi_2_Dec.begin(); itDec!=Pi_2_Dec.end(); ++itDec){
      Spin lamPi2=0;      
      boost::shared_ptr<const JPCLS> theDecState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));
      
      complex<double> tmpDec(0.,0.);
      
      tmpDec+=sqrt(2.*theDecState->L+1)*Clebsch(theDecState->L, 0, theDecState->S, 0, theDecState->J, 0)*Clebsch(0,0, 0, 0, theDecState->S, 0)
	*(BreitWigner(KKPi0, Pi_2_Mass, Pi_2_Width)*theData->DfPi2Tof2Pi0[theDecState->J][0][0]
	  + BreitWigner(KKPi1, Pi_2_Mass, Pi_2_Width)*theData->DfPi2Tof2Pi1[theDecState->J][0][0]);      

      currentResultDec+=tmpDec*theMagDec*expiphiDec;
    }

    result+=theMagProd*expiphiProd*currentResultDec;
  }

  result*=BreitWigner(KK, f0Mass, f0Width);
  return result;

}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK2KToK2PiAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_2_Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K_2_Dec, double K_2prod_Mass, double K_2prod_Width, double K_2decMass, double K_2decWidth){

  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
		     theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
		     theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
		     theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
		     theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
		     theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
		     theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itProd;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;

  for ( itProd=ChiToK_2_Prod.begin(); itProd!=ChiToK_2_Prod.end(); ++itProd){

    boost::shared_ptr<const JPCLS> theState=itProd->first;
    double theMagProd=itProd->second.first;
    double thePhiProd=itProd->second.second;
    complex<double> expiphiProd(cos(thePhiProd), sin(thePhiProd));

    complex<double> currentResultDec(0.,0.);

    Spin lamK2prod=0;

    for ( itDec=K_2_Dec.begin(); itDec!=K_2_Dec.end(); ++itDec){

      boost::shared_ptr<const JPCLS> theDecState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));
      complex<double> tmpDecKpPiPi(0.,0.);
      complex<double> tmpDecKmPiPi(0.,0.);

      for (Spin lamK2dec=-2; lamK2dec<=2; ++lamK2dec){
	tmpDecKpPiPi+=sqrt(2.*theDecState->L+1.)*Clebsch(theDecState->L, 0, theDecState->S, lamK2dec, theDecState->J, lamK2dec)
	  *Clebsch(2,lamK2dec, 0, 0, theDecState->S, lamK2dec)*sqrt(5.)
	  * BreitWigner(KpPiPi, K_2prod_Mass, K_2prod_Width)
	      *( BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K_2decMass, K_2decWidth,2)
		 *theData->DfKp2ToKppi1[theDecState->J][lamK2prod][lamK2dec]
		 *theData->DfK2pToKpPi0ViaKpPiPi[2][lamK2dec][0]
		 +
		 BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K_2decMass, K_2decWidth,2)
		 *theData->DfKp2ToKppi0[theDecState->J][lamK2prod][lamK2dec]
		 *theData->DfK2pToKpPi1ViaKpPiPi[2][lamK2dec][0] 
		 ); 

	tmpDecKmPiPi+=sqrt(2.*theDecState->L+1.)*Clebsch(theDecState->L, 0, theDecState->S, lamK2dec, theDecState->J, lamK2dec)
	  *Clebsch(2,lamK2dec, 0, 0, theDecState->S, lamK2dec)*sqrt(5.)
	  * BreitWigner(KmPiPi, K_2prod_Mass, K_2prod_Width)
	      *( BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K_2decMass, K_2decWidth,2)
		 *theData->DfKm2ToKmpi1[theDecState->J][lamK2prod][lamK2dec]
		 *theData->DfK2mToKmPi0ViaKmPiPi[2][lamK2dec][0]
		 +
		 BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K_2decMass, K_2decWidth,2)
		 *theData->DfKm2ToKmpi0[theDecState->J][lamK2prod][lamK2dec]
		 *theData->DfK2mToKmPi1ViaKmPiPi[2][lamK2dec][0] 
		 );


      }
      currentResultDec+=theMagDec*expiphiDec*(tmpDecKpPiPi+tmpDecKmPiPi);
    }

    result+=theMagProd*expiphiProd*currentResultDec;
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK1Tof980_piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1Dec, double K1Mass, double K1Width, double f980_Mass, double f980_gKK, double f980_gPiPi){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);



  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());

  Spin lamK1=0;
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itProd;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itDec;

  for ( itProd=K1Prod.begin(); itProd!=K1Prod.end(); ++itProd){

    boost::shared_ptr<const JPCLS> theState=itProd->first;
    double theMagProd=itProd->second.first;
    double thePhiProd=itProd->second.second;
    complex<double> expiphiProd(cos(thePhiProd), sin(thePhiProd));

    complex<double> currentResultDec(0.,0.);

    for ( itDec=K1Dec.begin(); itDec!=K1Dec.end(); ++itDec){
      
      boost::shared_ptr<const JPCLS> theDecState=itDec->first;
      double theMagDec=itDec->second.first;
      double thePhiDec=itDec->second.second;
      complex<double> expiphiDec(cos(thePhiDec), sin(thePhiDec));

      complex<double> ampKpPiPi=conj(theData->DfK1pTof0Kp[1][0][0])*BreitWigner(KpPiPi, K1Mass, K1Width);
      complex<double> ampKmPiPi=conj(theData->DfK1mTof0Km[1][0][0])*BreitWigner(KmPiPi, K1Mass, K1Width);


      currentResultDec+=theMagDec*expiphiDec*sqrt(2.*theState->L+1.)*Clebsch(theState->L, 0., theState->S, 0, theState->J, 0)*Clebsch(1,0, 0, 0, theState->S, 0)*Flatte(PiPi, decPair1, decPair2, f980_Mass, f980_gPiPi, f980_gKK)*(ampKpPiPi+ampKmPiPi);
    }

    result+=theMagProd*expiphiProd*currentResultDec;
  }
  return result;
}



void AbsPsi2STo2K2PiGamLh::print(std::ostream& os) const{
  os << "AbsPsi2STo2K2PiGamLh::print\n";
}

