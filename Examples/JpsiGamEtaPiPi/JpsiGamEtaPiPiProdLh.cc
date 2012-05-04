#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiProdLh.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiEventList.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiFitParams.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


JpsiGamEtaPiPiProdLh::JpsiGamEtaPiPiProdLh(boost::shared_ptr<const EvtDataBaseList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsLh(theEvtList)
  ,_etaToPiPiEtaHyp(false)
  ,_etaToa980PiHyp(false)
  ,_etaToa2_1320PiHyp(false)
  ,_etaToSigmaEtaHyp(false)
  ,_etaTof0_980EtaHyp(false)
  ,_etaTof2_1270EtaHyp(false)
  ,_eta2ToPiPiEtaHyp(false)
  ,_eta2Toa980PiHyp(false)
  ,_eta2Toa2_1320PiHyp(false)
  ,_eta2ToSigmaEtaHyp(false)
  ,_eta2Tof0_980EtaHyp(false)
  ,_eta2Tof2_1270EtaHyp(false)
  ,_f1ToPiPiEtaHyp(false)
  ,_f1Toa980PiHyp(false)
  ,_f1Toa2_1320PiHyp(false)
  ,_f1ToSigmaEtaHyp(false)
  ,_f1Tof0_980EtaHyp(false)
  ,_f1Tof2_1270EtaHyp(false)
  ,_usePhasespace(false)
  ,_massPi0(0.1349766)
  ,_massKplus(0.493677)
  ,_massK0(0.497614)
  ,_massEta(0.547853)
  ,_decPairPiPi(make_pair(_massPi0, _massPi0))
  ,_decPairKK(make_pair(_massKplus,_massK0))
  ,_decPairPi0Eta(make_pair(_massPi0, _massEta))
{
  initializeHypothesisMap( hypMap);
 
}

JpsiGamEtaPiPiProdLh::JpsiGamEtaPiPiProdLh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsLh(theLhPtr->getEventList())
  ,_etaToPiPiEtaHyp(false)
  ,_etaToa980PiHyp(false)
  ,_etaToa2_1320PiHyp(false)
  ,_etaToSigmaEtaHyp(false)
  ,_etaTof0_980EtaHyp(false)
  ,_etaTof2_1270EtaHyp(false)
  ,_eta2ToPiPiEtaHyp(false)
  ,_eta2Toa980PiHyp(false)
  ,_eta2Toa2_1320PiHyp(false)
  ,_eta2ToSigmaEtaHyp(false)
  ,_eta2Tof0_980EtaHyp(false)
  ,_eta2Tof2_1270EtaHyp(false)
  ,_f1ToPiPiEtaHyp(false)
  ,_f1Toa980PiHyp(false)
  ,_f1Toa2_1320PiHyp(false)
  ,_f1ToSigmaEtaHyp(false)
  ,_f1Tof0_980EtaHyp(false)
  ,_f1Tof2_1270EtaHyp(false)
  ,_usePhasespace(false)
  ,_massPi0(0.1349766)
  ,_massKplus(0.493677)
  ,_massK0(0.497614)
  ,_massEta(0.547853)
  ,_decPairPiPi(make_pair(_massPi0, _massPi0))
  ,_decPairKK(make_pair(_massKplus,_massK0))
  ,_decPairPi0Eta(make_pair(_massPi0, _massEta))
{
  
  initializeHypothesisMap( hypMap);
  
}

JpsiGamEtaPiPiProdLh::~JpsiGamEtaPiPiProdLh()
{;
}



double JpsiGamEtaPiPiProdLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;
  
 
  complex<double> JmpGmp(0.,0.);
  complex<double> JmpGmm(0.,0.);
  complex<double> JmmGmp(0.,0.);
  complex<double> JmmGmm(0.,0.);
  
  if(_etaToPiPiEtaHyp || _etaToa980PiHyp || _etaToa2_1320PiHyp ||_etaToSigmaEtaHyp || _etaToSigmaEtaHyp || _etaTof0_980EtaHyp || _etaTof2_1270EtaHyp){
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEtaGamMag=theParamVal.MagLamLams[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma];
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEtaGamPhi=theParamVal.PhiLamLams[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma];

    std::map<Spin,std::map<Spin, double > > MagProdMap;
    std::map<Spin,std::map<Spin, double > > PhiProdMap;
    std::map<Spin,std::map<Spin, double > > ParityProdMap;

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=PsiToEtaGamMag.begin(); itMag!=PsiToEtaGamMag.end(); ++itMag){
      boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
      MagProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=1.;
      PhiProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=PsiToEtaGamPhi[currentJPClamlam];

      MagProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=currentJPClamlam->parityFactor;
      PhiProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=PsiToEtaGamPhi[currentJPClamlam];      
    }

    complex<double>  JpGpTmp = ParityProdMap[0][1]*psiToXGammaAmp(1, 0, 0, 1, theData, MagProdMap[0][1], PhiProdMap[0][1]);
    complex<double>  JpGmTmp = ParityProdMap[0][-1]*psiToXGammaAmp(1, 0, 0, -1, theData, MagProdMap[0][-1], PhiProdMap[0][-1]);
    complex<double>  JmGpTmp = ParityProdMap[0][1]*psiToXGammaAmp(-1, 0, 0, 1, theData, MagProdMap[0][1], PhiProdMap[0][1]);
    complex<double>  JmGmTmp = ParityProdMap[0][-1]*psiToXGammaAmp(-1, 0, 0, -1, theData, MagProdMap[0][-1], PhiProdMap[0][-1]);
    complex<double>  TmpDecAmp(0.,0.);    

    if(_etaToPiPiEtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToPiPiEtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::EtaToPiPiEta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToPiPiEtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::EtaToPiPiEta];
      
	TmpDecAmp+=XToPiPiEtaAmp(0, 0, theData, etaToPiPiEtaMag, etaToPiPiEtaPhi);
    }

    if(_etaToa980PiHyp){

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToA980PiMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::EtaToA980Pi];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToA980PiPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::EtaToA980Pi];
      Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());

      Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());
  
      complex<double> a980PlusFlatte=A980DecFlatte(theParamVal, p4EtaPiplus);
      complex<double> a980MinusFlatte=A980DecFlatte(theParamVal, p4EtaPiminus);
      
	TmpDecAmp+=XToAPiAmp(0, 0, 0, theData, etaToA980PiMag, etaToA980PiPhi, a980PlusFlatte, a980MinusFlatte);


    }

    if(_etaToa2_1320PiHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToA2_1320PiMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToA2_1320PiPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi];

      double a2_1320Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::a2_1320];
      double a2_1320Width=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::a2_1320];

      Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());
      
      Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				    theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());
 
     complex<double> a2_1320PlusBW=BreitWigner(p4EtaPiplus, a2_1320Mass, a2_1320Width);
     complex<double> a2_1320MinusBW=BreitWigner(p4EtaPiminus, a2_1320Mass, a2_1320Width);

      TmpDecAmp+=XToAPiAmp(0, 0, 2, theData, etaToA2_1320PiMag, etaToA2_1320PiPhi, a2_1320PlusBW, a2_1320MinusBW);
    }

    if(_etaToSigmaEtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToSigmaEtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::EtaToSigmaEta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaToSigmaEtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::EtaToSigmaEta];

      double sigmaMass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::sigma];
      double sigmaWidth=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::sigma];
      //dummy

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());
  
      complex<double> dynamicFct = BreitWigner(p4PiplusPiminus, sigmaMass, sigmaWidth);
 
     TmpDecAmp+=XToEtaFAmp(0, 0, 0, theData, etaToSigmaEtaMag, etaToSigmaEtaPhi, dynamicFct);
    }



    if(_etaTof0_980EtaHyp){

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaTof0_980EtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::EtaTof0_980Eta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaTof0_980EtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::EtaTof0_980Eta];

      double f0_980Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::f0_980];
      double f0_980gPiPi=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gPiPi];
      double f0_980gKK=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gKK];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());

      
      complex<double> dynamicFct=Flatte(p4PiplusPiminus, _decPairPiPi, _decPairKK, f0_980Mass, f0_980gPiPi, f0_980gKK);

      TmpDecAmp+=XToEtaFAmp(0, 0, 0, theData, etaTof0_980EtaMag, etaTof0_980EtaPhi, dynamicFct);

    } 

   if(_etaTof2_1270EtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaTof2_1270EtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::EtaTof2_1270Eta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > etaTof2_1270EtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::EtaTof2_1270Eta];

      double f2_1270Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::f2_1270];
      double f2_1270Width=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::f2_1270];
      //dummy

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());
  
      complex<double> dynamicFct = BreitWigner(p4PiplusPiminus, f2_1270Mass, f2_1270Width);
 
      TmpDecAmp+=XToEtaFAmp(0, 0, 2, theData, etaTof2_1270EtaMag, etaTof2_1270EtaPhi, dynamicFct);
    }

      JmpGmp+=JpGpTmp*TmpDecAmp;
      JmpGmm+=JpGmTmp*TmpDecAmp;
      JmmGmp+=JmGpTmp*TmpDecAmp;
      JmmGmm+=JmGmTmp*TmpDecAmp;


    
  }

  
  if(_f1ToPiPiEtaHyp || _f1Toa980PiHyp || _f1Toa2_1320PiHyp || _f1ToSigmaEtaHyp || _f1Tof0_980EtaHyp || _f1Tof2_1270EtaHyp){
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiTof1GamMag=theParamVal.MagLamLams[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma];
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiTof1GamPhi=theParamVal.PhiLamLams[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma];

    std::map<Spin,std::map<Spin, double > > MagProdMap;
    std::map<Spin,std::map<Spin, double > > PhiProdMap;
    std::map<Spin,std::map<Spin, double > > ParityProdMap;

    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;
    for (itMag=PsiTof1GamMag.begin(); itMag!=PsiTof1GamMag.end(); ++itMag){
      boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
      MagProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=1.;
      PhiProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=PsiTof1GamPhi[currentJPClamlam];

      MagProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=currentJPClamlam->parityFactor;
      PhiProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=PsiTof1GamPhi[currentJPClamlam];      
    }

    std::map<Spin,complex<double> > JpGpTmpMap;
    std::map<Spin,complex<double> > JpGmTmpMap;
    std::map<Spin,complex<double> > JmGpTmpMap;
    std::map<Spin,complex<double> > JmGmTmpMap;
    std::map<Spin,complex<double> > TmpDecAmp;
 
    for (Spin helf1=0; helf1<2; helf1++){
      JpGpTmpMap[helf1]= ParityProdMap[helf1][1]*psiToXGammaAmp(1, 1, helf1, 1, theData, MagProdMap[helf1][1], PhiProdMap[helf1][1]);
      JpGmTmpMap[-helf1]= ParityProdMap[-helf1][-1]*psiToXGammaAmp(1, 1, -helf1, -1, theData, MagProdMap[-helf1][-1], PhiProdMap[-helf1][-1]); 
      JmGpTmpMap[helf1]= ParityProdMap[helf1][1]*psiToXGammaAmp(-1, 1, helf1, 1, theData, MagProdMap[helf1][1], PhiProdMap[helf1][1]);
      JmGmTmpMap[-helf1]= ParityProdMap[-helf1][-1]*psiToXGammaAmp(-1, 1, -helf1, -1, theData, MagProdMap[-helf1][-1], PhiProdMap[-helf1][-1]);
      TmpDecAmp[helf1] = complex<double> (0.,0.);
    }
  

    if(_f1ToPiPiEtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToPiPiEtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::F1ToPiPiEta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToPiPiEtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::F1ToPiPiEta];
      
      for (Spin helf1=-1; helf1<2; helf1++){
	TmpDecAmp[helf1]+=XToPiPiEtaAmp(1, helf1, theData, f1ToPiPiEtaMag, f1ToPiPiEtaPhi);
      }
      
    }

    if(_f1Toa980PiHyp){

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToA980PiMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::F1ToA980Pi];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToA980PiPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::F1ToA980Pi];
      Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());

      Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());
  
      complex<double> a980PlusFlatte=A980DecFlatte(theParamVal, p4EtaPiplus);
      complex<double> a980MinusFlatte=A980DecFlatte(theParamVal, p4EtaPiminus);
      
      for (Spin helf1=-1; helf1<2; helf1++){
	TmpDecAmp[helf1]+=XToAPiAmp(1, helf1, 0, theData, f1ToA980PiMag, f1ToA980PiPhi, a980PlusFlatte, a980MinusFlatte);
      }

    }



    if(_f1Toa2_1320PiHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToA2_1320PiMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::F1ToA2_1320Pi];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToA2_1320PiPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::F1ToA2_1320Pi];

      double a2_1320Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::a2_1320];
      double a2_1320Width=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::a2_1320];

      Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());
      
      Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				    theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());
 
     complex<double> a2_1320PlusBW=BreitWigner(p4EtaPiplus, a2_1320Mass, a2_1320Width);
     complex<double> a2_1320MinusBW=BreitWigner(p4EtaPiminus, a2_1320Mass, a2_1320Width);

     for (Spin helif1=-1; helif1<=1; helif1++){
       TmpDecAmp[helif1]+=XToAPiAmp(1, helif1, 2, theData, f1ToA2_1320PiMag, f1ToA2_1320PiPhi, a2_1320PlusBW, a2_1320MinusBW);
     }
    }



    if(_f1ToSigmaEtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToSigmaEtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::F1ToSigmaEta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1ToSigmaEtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::F1ToSigmaEta];

      double sigmaMass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::sigma];
      double sigmaWidth=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::sigma];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());
  
      complex<double> dynamicFct = BreitWigner(p4PiplusPiminus, sigmaMass, sigmaWidth);

      for (Spin helif1=-1; helif1<=1; helif1++){
	TmpDecAmp[helif1] += XToEtaFAmp( 1, helif1, 0, theData, f1ToSigmaEtaMag, f1ToSigmaEtaPhi, dynamicFct);
      }
    }

    if(_f1Tof0_980EtaHyp){

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1Tof0_980EtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::F1Tof0_980Eta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1Tof0_980EtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::F1Tof0_980Eta];

      double f0_980Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::f0_980];
      double f0_980gPiPi=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gPiPi];
      double f0_980gKK=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gKK];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());

      
      complex<double> dynamicFct=Flatte(p4PiplusPiminus, _decPairPiPi, _decPairKK, f0_980Mass, f0_980gPiPi, f0_980gKK);

      for (Spin helif1=-1; helif1<=1; helif1++){ 
	TmpDecAmp[helif1] += XToEtaFAmp(1, helif1, 0, theData, f1Tof0_980EtaMag, f1Tof0_980EtaPhi, dynamicFct);
      }

    }

   if(_f1Tof2_1270EtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1Tof2_1270EtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::F1Tof2_1270Eta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > f1Tof2_1270EtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::F1Tof2_1270Eta];

      double f2_1270Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::f2_1270];
      double f2_1270Width=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::f2_1270];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());
  
      complex<double> dynamicFct = BreitWigner(p4PiplusPiminus, f2_1270Mass, f2_1270Width);

      for (Spin helif1=-1; helif1<=1; helif1++){
	TmpDecAmp[helif1] += XToEtaFAmp( 1, helif1, 2, theData, f1Tof2_1270EtaMag, f1Tof2_1270EtaPhi, dynamicFct);
      }
    }



    for (Spin helf1=-1; helf1<2; helf1++){
      JmpGmp+=JpGpTmpMap[helf1]*TmpDecAmp[helf1];
      JmpGmm+=JpGmTmpMap[helf1]*TmpDecAmp[helf1];
      JmmGmp+=JmGpTmpMap[helf1]*TmpDecAmp[helf1];
      JmmGmm+=JmGmTmpMap[helf1]*TmpDecAmp[helf1];

    } 
  }







  if(_eta2ToPiPiEtaHyp || _eta2Toa980PiHyp || _eta2Toa2_1320PiHyp ||_etaToSigmaEtaHyp || _etaTof0_980EtaHyp || _etaTof2_1270EtaHyp){
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEta2GamMag=theParamVal.MagLamLams[paramEnumJpsiGamEtaPiPi::PsiToEta2Gamma];
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > PsiToEta2GamPhi=theParamVal.PhiLamLams[paramEnumJpsiGamEtaPiPi::PsiToEta2Gamma];

    std::map<Spin,std::map<Spin, double > > MagProdMap;
    std::map<Spin,std::map<Spin, double > > PhiProdMap;
    std::map<Spin,std::map<Spin, double > > ParityProdMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess >::iterator itMag;

    for (itMag=PsiToEta2GamMag.begin(); itMag!=PsiToEta2GamMag.end(); ++itMag){
      boost::shared_ptr<const JPClamlam> currentJPClamlam=itMag->first;
      MagProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=1.;
      PhiProdMap[currentJPClamlam->lam1][currentJPClamlam->lam2]=PsiToEta2GamPhi[currentJPClamlam];

      MagProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=itMag->second;
      ParityProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=currentJPClamlam->parityFactor;
      PhiProdMap[-currentJPClamlam->lam1][-currentJPClamlam->lam2]=PsiToEta2GamPhi[currentJPClamlam];      
    }

    std::map<Spin,complex<double> > JpGpTmpMap;
    std::map<Spin,complex<double> > JpGmTmpMap;
    std::map<Spin,complex<double> > JmGpTmpMap;
    std::map<Spin,complex<double> > JmGmTmpMap;
    std::map<Spin,complex<double> > TmpDecAmp;
 
    for (Spin heliEta2=0; heliEta2<=2; heliEta2++){
      JpGpTmpMap[heliEta2]= ParityProdMap[heliEta2][1]*psiToXGammaAmp(1, 2, heliEta2, 1, theData, MagProdMap[heliEta2][1], PhiProdMap[heliEta2][1]);
      JpGmTmpMap[-heliEta2]= ParityProdMap[-heliEta2][-1]*psiToXGammaAmp(1, 2, -heliEta2, -1, theData, MagProdMap[-heliEta2][-1], PhiProdMap[-heliEta2][-1]); 
      JmGpTmpMap[heliEta2]=ParityProdMap[heliEta2][1]*psiToXGammaAmp(-1, 2, heliEta2, 1, theData, MagProdMap[heliEta2][1], PhiProdMap[heliEta2][1]);
      JmGmTmpMap[-heliEta2]=ParityProdMap[-heliEta2][-1]*psiToXGammaAmp(-1, 2, -heliEta2, -1, theData, MagProdMap[-heliEta2][-1], PhiProdMap[-heliEta2][-1]);
      TmpDecAmp[heliEta2] = complex<double> (0.,0.);
    }
  

    if(_eta2ToPiPiEtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToPiPiEtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::Eta2ToPiPiEta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToPiPiEtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::Eta2ToPiPiEta];
      
      for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){
	TmpDecAmp[heliEta2]+=XToPiPiEtaAmp(2, heliEta2, theData, eta2ToPiPiEtaMag, eta2ToPiPiEtaPhi);
      }
      
    }

    if(_eta2Toa980PiHyp){

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToA980PiMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::Eta2ToA980Pi];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToA980PiPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::Eta2ToA980Pi];
      Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());

      Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());
  
      complex<double> a980PlusFlatte=A980DecFlatte(theParamVal, p4EtaPiplus);
      complex<double> a980MinusFlatte=A980DecFlatte(theParamVal, p4EtaPiminus);
      
      for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){
	TmpDecAmp[heliEta2]+=XToAPiAmp(2, heliEta2, 0, theData, eta2ToA980PiMag, eta2ToA980PiPhi, a980PlusFlatte, a980MinusFlatte);
      }

    }


    if(_eta2Toa2_1320PiHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToA2_1320PiMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::Eta2ToA2_1320Pi];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToA2_1320PiPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::Eta2ToA2_1320Pi];

      double a2_1320Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::a2_1320];
      double a2_1320Width=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::a2_1320];

      Vector4<double > p4EtaPiplus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].E(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPip_HeliPsi].Pz());
      
      Vector4<double > p4EtaPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].E(),
				    theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Px(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Py(),
				   theData->FourVecs[enumJpsiGamEtaPiPiData::V4_EtaPim_HeliPsi].Pz());
 
     complex<double> a2_1320PlusBW=BreitWigner(p4EtaPiplus, a2_1320Mass, a2_1320Width);
     complex<double> a2_1320MinusBW=BreitWigner(p4EtaPiminus, a2_1320Mass, a2_1320Width);

     for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){
       TmpDecAmp[heliEta2]+=XToAPiAmp(2, heliEta2, 2, theData, eta2ToA2_1320PiMag, eta2ToA2_1320PiPhi, a2_1320PlusBW, a2_1320MinusBW);
     }
    }


    if(_eta2ToSigmaEtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToSigmaEtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::Eta2ToSigmaEta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2ToSigmaEtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::Eta2ToSigmaEta];

      double sigmaMass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::sigma];
      double sigmaWidth=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::sigma];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());
  
      complex<double> dynamicFct = BreitWigner(p4PiplusPiminus, sigmaMass, sigmaWidth);

      for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){
	TmpDecAmp[heliEta2] += XToEtaFAmp( 2, heliEta2, 0, theData, eta2ToSigmaEtaMag, eta2ToSigmaEtaPhi, dynamicFct);
      }
    }

    if(_eta2Tof0_980EtaHyp){

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2Tof0_980EtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::Eta2Tof0_980Eta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2Tof0_980EtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::Eta2Tof0_980Eta];

      double f0_980Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::f0_980];
      double f0_980gPiPi=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gPiPi];
      double f0_980gKK=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gKK];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());

      
      complex<double> dynamicFct=Flatte(p4PiplusPiminus, _decPairPiPi, _decPairKK, f0_980Mass, f0_980gPiPi, f0_980gKK);

      for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){ 
	TmpDecAmp[heliEta2] += XToEtaFAmp(2, heliEta2, 0, theData, eta2Tof0_980EtaMag, eta2Tof0_980EtaPhi, dynamicFct);
      }

    }
 
   if(_eta2Tof2_1270EtaHyp){
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2Tof2_1270EtaMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::Eta2Tof2_1270Eta];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > eta2Tof2_1270EtaPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::Eta2Tof2_1270Eta];

      double f2_1270Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::f2_1270];
      double f2_1270Width=theParamVal.Widths[paramEnumJpsiGamEtaPiPi::f2_1270];

      Vector4<double > p4PiplusPiminus(theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].E(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Px(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Py(),
				       theData->FourVecs[enumJpsiGamEtaPiPiData::V4_PipPim_HeliPsi].Pz());
  
      complex<double> dynamicFct = BreitWigner(p4PiplusPiminus, f2_1270Mass, f2_1270Width);

      for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){
	TmpDecAmp[heliEta2] += XToEtaFAmp( 2, heliEta2, 2, theData, eta2Tof2_1270EtaMag, eta2Tof2_1270EtaPhi, dynamicFct);
      }
    }



    for (Spin heliEta2=-2; heliEta2<=2; heliEta2++){
      JmpGmp+=JpGpTmpMap[heliEta2]*TmpDecAmp[heliEta2];
      JmpGmm+=JpGmTmpMap[heliEta2]*TmpDecAmp[heliEta2];
      JmmGmp+=JmGpTmpMap[heliEta2]*TmpDecAmp[heliEta2];
      JmmGmm+=JmGmTmpMap[heliEta2]*TmpDecAmp[heliEta2];

    } 
  }






  

  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);

  if(_usePhasespace){
    result = result + theParamVal.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace];
  }
  
  return result;  
}

complex<double> JpsiGamEtaPiPiProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}

complex<double> JpsiGamEtaPiPiProdLh::psiToXGammaAmp(Spin Minit, Spin jX, Spin lamX, Spin lamGamma, EvtData* theData, 
						     double PsiToXGamMag, double PsiToXGamPhi ){
   Spin lambda = lamX-lamGamma;
   complex<double> expiphiPsi(cos(PsiToXGamPhi), sin(PsiToXGamPhi));
   complex<double> result = PsiToXGamMag*expiphiPsi*conj( theData->WignerDs[enumJpsiGamEtaPiPiData::Df_Psi][1][Minit][lambda]  );

   return result;

}

complex<double> JpsiGamEtaPiPiProdLh::XToPiPiEtaAmp(Spin jX, Spin lamX, EvtData* theData, 
						    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaMag, 
						    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToPiPiEtaPhi){
   complex<double> result(0.,0.);

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;
   for ( itXMag=XToPiPiEtaMag.begin(); itXMag!=XToPiPiEtaMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToPiPiEtaPhi[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));

     complex<double> amp = theXMag*expiphiX*sqrt(2.*XState->L+1.)
       *Clebsch(XState->L, 0, XState->S, 0, XState->J, 0)
       *Clebsch(0, 0, 0, 0, XState->S, 0);
       
     result+= amp;
   }
   result*=conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_etapipidec][jX][lamX][0]);
   return result;
}

complex<double> JpsiGamEtaPiPiProdLh::XToAPiAmp(Spin jX, Spin lamX, Spin jA, EvtData* theData, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiMag, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToAPiPhi,
						complex<double>& dynAplus, complex<double>& dynAminus){

   complex<double> result(0.,0.);

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

   for ( itXMag=XToAPiMag.begin(); itXMag!=XToAPiMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToAPiPhi[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
     complex<double> amp(0.,0.);     
     for(Spin lamA = -jA; lamA <= jA; lamA++){
       
       if( fabs(lamA)> XState->J || fabs(lamA)>XState->S) continue;
       
       amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)
	 *Clebsch(XState->L, 0, XState->S, lamA, XState->J, lamA)
	 *Clebsch(jA, lamA, 0, 0, XState->S, lamA)*
	 ( conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_XToAplusPiminusdec][jX][lamX][lamA])
	   *conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_AplusDec][jA][lamA][0])*dynAplus
	  + conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_XToAminusPiplusdec][jX][lamX][lamA])
	   *conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_AminusDec][jA][lamA][0])*dynAminus );
     }
     
     result+= amp;
   }
   
   return result;
}


complex<double> JpsiGamEtaPiPiProdLh::XToEtaFAmp(Spin jX, Spin lamX, Spin jf, EvtData* theData, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFMag, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFPhi, 
						 complex<double> dynamicFct){
  complex<double> result(0.,0.);
  
  
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itXMag;

   for ( itXMag=XToEtaFMag.begin(); itXMag!=XToEtaFMag.end(); ++itXMag){
     boost::shared_ptr<const JPCLS> XState=itXMag->first;
     double theXMag=itXMag->second;
     double theXPhi=XToEtaFPhi[XState];
     complex<double> expiphiX(cos(theXPhi), sin(theXPhi));
     complex<double> amp(0.,0.);     
     for(Spin lamf = -jf; lamf <= jf; lamf++){
       
       if( fabs(lamf)> XState->J || fabs(lamf)>XState->S) continue;
       
       amp += theXMag*expiphiX*sqrt(2.*XState->L+1.)
	 *Clebsch(XState->L, 0, XState->S, lamf, XState->J, lamf)
	 *Clebsch(jf, lamf, 0, 0, XState->S, lamf)
	 *conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_XTofEtadec][jX][lamX][lamf])
	 *conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_fToPiPiDec][jf][lamf][0])*dynamicFct;
     }
     
     result+= amp;
   }
   
   return result;
}




complex<double> JpsiGamEtaPiPiProdLh::A980DecFlatte(fitParams& theParamVal, const Vector4<double> &__p4){
  
  complex<double> result(0.,0.);


  double a0_980Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::a0_980];
  double a0_980gPiEta=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gPiEta];
  double a0_980gKK=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gKK];

  result+=Flatte(__p4, _decPairPi0Eta, _decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK);

  return result;
}


void JpsiGamEtaPiPiProdLh::print(std::ostream& os) const{
  os << "JpsiGamEtaPiPiProdLh::print\n";
}


void JpsiGamEtaPiPiProdLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  JpsiGamEtaPiPiFitParams theFitParams;

  std::map<int, std::vector< boost::shared_ptr<const JPClamlam> > > theAmpLamLamMap;  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;
  
  if(_etaToPiPiEtaHyp || _etaToa980PiHyp || _etaToa2_1320PiHyp || _etaToSigmaEtaHyp || _etaTof0_980EtaHyp || _etaTof2_1270EtaHyp){

    theAmpLamLamMap[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma] = theFitParams.jpcLamLamVec(paramEnumJpsiGamEtaPiPi::PsiToEtaGamma);

    if(_etaToPiPiEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToPiPiEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToPiPiEta);
    }

    if (_etaToa980PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToA980Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToA980Pi);
    }
    
    if(_etaToa2_1320PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi);
    }

    if(_etaToSigmaEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToSigmaEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToSigmaEta);
    }

    if(_etaTof0_980EtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaTof0_980Eta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaTof0_980Eta);
    }

    if(_etaTof2_1270EtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaTof2_1270Eta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaTof2_1270Eta);
    }
  }


 if(_eta2ToPiPiEtaHyp || _eta2Toa980PiHyp || _eta2Toa2_1320PiHyp || _eta2ToSigmaEtaHyp || _eta2Tof0_980EtaHyp || _eta2Tof2_1270EtaHyp){
    theAmpLamLamMap[paramEnumJpsiGamEtaPiPi::PsiToEta2Gamma] = theFitParams.jpcLamLamVec(paramEnumJpsiGamEtaPiPi::PsiToEta2Gamma);

    if(_eta2ToPiPiEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::Eta2ToPiPiEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::Eta2ToPiPiEta);
    }

    if (_eta2Toa980PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::Eta2ToA980Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::Eta2ToA980Pi);
    }
    
    if(_eta2Toa2_1320PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::Eta2ToA2_1320Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::Eta2ToA2_1320Pi);
    }

    if(_eta2ToSigmaEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::Eta2ToSigmaEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::Eta2ToSigmaEta);
    }

    if(_eta2Tof0_980EtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::Eta2Tof0_980Eta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::Eta2Tof0_980Eta);
    }

    if(_eta2Tof2_1270EtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::Eta2Tof2_1270Eta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::Eta2Tof2_1270Eta);
    }
  }



  
  if(_f1ToPiPiEtaHyp || _f1Toa980PiHyp|| _f1Toa2_1320PiHyp || _f1ToSigmaEtaHyp || _f1Tof0_980EtaHyp || _f1Tof2_1270EtaHyp){
    theAmpLamLamMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma] = theFitParams.jpcLamLamVec(paramEnumJpsiGamEtaPiPi::PsiToF1Gamma);
    
    if(_f1ToPiPiEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1ToPiPiEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1ToPiPiEta);
    }

    if (_f1Toa980PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1ToA980Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1ToA980Pi);
    }

    if(_f1Toa2_1320PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1ToA2_1320Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1ToA2_1320Pi);
    }

    if(_f1ToSigmaEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1ToSigmaEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1ToSigmaEta);
    }

    if(_f1Tof0_980EtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1Tof0_980Eta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1Tof0_980Eta);
    }

    if(_f1Tof2_1270EtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1Tof2_1270Eta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1Tof2_1270Eta);
    }
  }
  
  std::map<int, std::vector< boost::shared_ptr<const JPClamlam> > >::iterator itAmpLamLamMap;
  for (itAmpLamLamMap=theAmpLamLamMap.begin(); itAmpLamLamMap!=theAmpLamLamMap.end(); ++itAmpLamLamMap){
    
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > valMagMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > errMagMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > valPhiMap;
    std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > errPhiMap;
    
    std::vector< boost::shared_ptr<const JPClamlam> >::iterator itAmp;
    for (itAmp=itAmpLamLamMap->second.begin(); itAmp!=itAmpLamLamMap->second.end(); ++itAmp){
      valMagMap[(*itAmp)]=0.1;
      errMagMap[(*itAmp)]=0.9; 
      valPhiMap[(*itAmp)]=0.0;
      errPhiMap[(*itAmp)]=0.8;      
    }
    
    fitVal.MagLamLams[itAmpLamLamMap->first]=valMagMap;
    fitVal.PhiLamLams[itAmpLamLamMap->first]=valPhiMap;  
    fitErr.MagLamLams[itAmpLamLamMap->first]=errMagMap;
    fitErr.PhiLamLams[itAmpLamLamMap->first]=errPhiMap;  
  }
  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > >::iterator itAmpMap;
  for (itAmpMap=theAmpMap.begin(); itAmpMap!=theAmpMap.end(); ++itAmpMap){
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > valMagMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > errMagMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > valPhiMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > errPhiMap;
    
    std::vector< boost::shared_ptr<const JPCLS> >::iterator itAmp;
    for (itAmp=itAmpMap->second.begin(); itAmp!=itAmpMap->second.end(); ++itAmp){
      valMagMap[(*itAmp)]=0.1;
      errMagMap[(*itAmp)]=0.9; 
      valPhiMap[(*itAmp)]=0.0;
      errPhiMap[(*itAmp)]=0.8;      
    }
    
    fitVal.Mags[itAmpMap->first]=valMagMap;
    fitVal.Phis[itAmpMap->first]=valPhiMap;  
    fitErr.Mags[itAmpMap->first]=errMagMap;
    fitErr.Phis[itAmpMap->first]=errPhiMap;  
  }

  
  //fill masses and widths
  if (_f1Toa980PiHyp || _etaToa980PiHyp ){
    fitVal.Masses[paramEnumJpsiGamEtaPiPi::a0_980]=0.98;
    fitErr.Masses[paramEnumJpsiGamEtaPiPi::a0_980]=0.03;
    fitVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gPiEta]= 1.0;
    fitErr.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gPiEta]=0.06;
    fitVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gKK]= 1.0;
    fitErr.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gKK]=0.06;
  }

  if(_etaTof0_980EtaHyp || _eta2Tof0_980EtaHyp || _f1Tof0_980EtaHyp){
    fitVal.Masses[paramEnumJpsiGamEtaPiPi::f0_980]=0.98;
    fitErr.Masses[paramEnumJpsiGamEtaPiPi::f0_980]=0.03;
    fitVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gPiPi] = 0.7;
    fitErr.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gPiPi]= 0.06;
    fitVal.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gKK]= 1.2;
    fitErr.gFactors[paramEnumJpsiGamEtaPiPi::f0_980gKK]=0.06;
 }

  if(_etaToSigmaEtaHyp || _eta2ToSigmaEtaHyp || _f1ToSigmaEtaHyp){
    fitVal.Masses[paramEnumJpsiGamEtaPiPi::sigma]=0.5;
    fitErr.Masses[paramEnumJpsiGamEtaPiPi::sigma]=0.05;
    fitVal.Widths[paramEnumJpsiGamEtaPiPi::sigma]=0.4;
    fitErr.Widths[paramEnumJpsiGamEtaPiPi::sigma]=0.08;
 }

  if(_etaToa2_1320PiHyp || _eta2Toa2_1320PiHyp || _f1Toa2_1320PiHyp){
    fitVal.Masses[paramEnumJpsiGamEtaPiPi::a2_1320]=1.32;
    fitErr.Masses[paramEnumJpsiGamEtaPiPi::a2_1320]=0.02;
    fitVal.Widths[paramEnumJpsiGamEtaPiPi::a2_1320]=0.1;
    fitErr.Widths[paramEnumJpsiGamEtaPiPi::a2_1320]=0.03;
  }

  if(_etaTof2_1270EtaHyp || _eta2Tof2_1270EtaHyp || _f1Tof2_1270EtaHyp){
    fitVal.Masses[paramEnumJpsiGamEtaPiPi::f2_1270]=1.275;
    fitErr.Masses[paramEnumJpsiGamEtaPiPi::f2_1270]=0.03;
    fitVal.Widths[paramEnumJpsiGamEtaPiPi::f2_1270]=0.185;
    fitErr.Widths[paramEnumJpsiGamEtaPiPi::f2_1270]=0.06;
  }

  //fill other params  
  if(_usePhasespace){
    fitVal.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace]=0.2;
    fitErr.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace]=0.4;
  }

}
  



bool 
JpsiGamEtaPiPiProdLh::initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   ){


  std::string theKey="etaToPiPiEtaHyp";
  setHyps( hypMap, _etaToPiPiEtaHyp, theKey);
  
  theKey="etaToa980PiHyp";
  setHyps( hypMap, _etaToa980PiHyp, theKey);

  theKey="etaToa2_1320PiHyp";
  setHyps( hypMap, _etaToa2_1320PiHyp, theKey);

  theKey="etaToSigmaEtaHyp";
  setHyps( hypMap, _etaToSigmaEtaHyp, theKey);

  theKey="etaTof0_980EtaHyp";
  setHyps( hypMap, _etaTof0_980EtaHyp, theKey);

  theKey="etaTof2_1270EtaHyp";
  setHyps( hypMap, _etaTof2_1270EtaHyp, theKey);

  theKey="eta2ToPiPiEtaHyp";
  setHyps( hypMap, _eta2ToPiPiEtaHyp, theKey);

  theKey="eta2Toa980PiHyp";
  setHyps( hypMap, _eta2Toa980PiHyp, theKey);

  theKey="eta2Toa2_1320PiHyp";
  setHyps( hypMap, _eta2Toa2_1320PiHyp, theKey);

  theKey="eta2ToSigmaEtaHyp";
  setHyps( hypMap, _eta2ToSigmaEtaHyp, theKey);

  theKey="eta2Tof0_980EtaHyp";
  setHyps( hypMap, _eta2Tof0_980EtaHyp, theKey);

  theKey="eta2Tof2_1270EtaHyp";
  setHyps( hypMap, _eta2Tof2_1270EtaHyp, theKey);

  theKey="f1ToPiPiEtaHyp";
  setHyps( hypMap, _f1ToPiPiEtaHyp, theKey);

  theKey="f1Toa980PiHyp";
  setHyps( hypMap, _f1Toa980PiHyp, theKey);

  theKey="f1Toa2_1320PiHyp";
  setHyps( hypMap, _f1Toa2_1320PiHyp, theKey);

  theKey="f1ToSigmaEtaHyp";
  setHyps( hypMap, _f1ToSigmaEtaHyp, theKey);

  theKey="f1Tof0_980EtaHyp";
  setHyps( hypMap, _f1Tof0_980EtaHyp, theKey);

  theKey="f1Tof2_1270EtaHyp";
  setHyps( hypMap, _f1Tof2_1270EtaHyp, theKey);

  theKey="usePhasespace";
  setHyps( hypMap, _usePhasespace, theKey);

  
  return true;
}








