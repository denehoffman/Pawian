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
  ,_f1ToPiPiEtaHyp(false)
  ,_f1Toa980PiHyp(false)
  ,_usePhasespace(false)
{
  initializeHypothesisMap( hypMap);
 
}

JpsiGamEtaPiPiProdLh::JpsiGamEtaPiPiProdLh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsLh(theLhPtr->getEventList())
  ,_etaToPiPiEtaHyp(false)
  ,_etaToa980PiHyp(false)
  ,_etaToa2_1320PiHyp(false)
  ,_f1ToPiPiEtaHyp(false)
  ,_f1Toa980PiHyp(false)
  ,_usePhasespace(false)
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
  
  if(_etaToPiPiEtaHyp || _etaToa980PiHyp || _etaToa2_1320PiHyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma];

    complex<double>  JpGpTmp = psiToXGammaAmp(1, 0, 0, 1, theData, PsiToEtacGamMag, PsiToEtacGamPhi);
    complex<double>  JpGmTmp = psiToXGammaAmp(1, 0, 0, -1, theData, PsiToEtacGamMag, PsiToEtacGamPhi);
    complex<double>  JmGpTmp = psiToXGammaAmp(-1, 0, 0, 1, theData, PsiToEtacGamMag, PsiToEtacGamPhi);
    complex<double>  JmGmTmp = psiToXGammaAmp(-1, 0, 0, -1, theData, PsiToEtacGamMag, PsiToEtacGamPhi);
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

      JmpGmp+=JpGpTmp*TmpDecAmp;
      JmpGmm+=JpGmTmp*TmpDecAmp;
      JmmGmp+=JmGpTmp*TmpDecAmp;
      JmmGmm+=JmGmTmp*TmpDecAmp;


    //    calcEtaGammaAmp( theData, PsiToEtacGamMag, PsiToEtacGamPhi, JmpGmp, JmpGmm, JmmGmp, JmmGmm);
    
  }

//   std::map<int, map<Spin,map<Spin,map<Spin,complex<double> > > > > f1Amps;
  if(_f1ToPiPiEtaHyp || _f1Toa980PiHyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof1GamMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof1GamPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma];

    std::map<Spin,complex<double> > JpGpTmpMap;
    std::map<Spin,complex<double> > JpGmTmpMap;
    std::map<Spin,complex<double> > JmGpTmpMap;
    std::map<Spin,complex<double> > JmGmTmpMap;
    std::map<Spin,complex<double> > TmpDecAmp;
 
    for (Spin helf1=-1; helf1<2; helf1++){
      JpGpTmpMap[helf1]= psiToXGammaAmp(1, 1, helf1, 1, theData, PsiTof1GamMag, PsiTof1GamPhi);
      JpGmTmpMap[helf1]= psiToXGammaAmp(1, 1, helf1, -1, theData, PsiTof1GamMag, PsiTof1GamPhi); 
      JmGpTmpMap[helf1]=psiToXGammaAmp(-1, 1, helf1, 1, theData, PsiTof1GamMag, PsiTof1GamPhi);
      JmGmTmpMap[helf1]=psiToXGammaAmp(-1, 1, helf1, -1, theData, PsiTof1GamMag, PsiTof1GamPhi);
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

    for (Spin helf1=-1; helf1<2; helf1++){
      JmpGmp+=JpGpTmpMap[helf1]*TmpDecAmp[helf1];
      JmpGmm+=JpGmTmpMap[helf1]*TmpDecAmp[helf1];
      JmmGmp+=JmGpTmpMap[helf1]*TmpDecAmp[helf1];
      JmmGmm+=JmGmTmpMap[helf1]*TmpDecAmp[helf1];

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
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToXGamMag, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToXGamPhi ){
   complex<double> result(0.,0.);
   Spin lambda = lamX-lamGamma;

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToXGamMag.begin(); itPsi!=PsiToXGamMag.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second;
     double thePsiPhi=PsiToXGamPhi[PsiState];
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));

     if( fabs(lambda)> PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(jX, lamX, 1, -lamGamma, PsiState->S, lambda  );
     
     result+= amp;
   }

   result*= conj( theData->WignerDs[enumJpsiGamEtaPiPiData::Df_Psi][1][Minit][lambda]  );
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
	 (conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_AplusPiminusdec][jX][lamX][lamA])*dynAplus+conj(theData->WignerDs[enumJpsiGamEtaPiPiData::Df_AminusPiplusdec][jX][lamX][lamA])*dynAminus);
     }
     
     result+= amp;
   }
   
   return result;
}


void JpsiGamEtaPiPiProdLh::calcEtaGammaAmp( EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaGamMag, 
						       std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaGamPhi, 
						       complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm){

  JmpGmp+=psiToXGammaAmp(1, 0, 0, 1, theData, PsiToEtaGamMag, PsiToEtaGamPhi );
  JmpGmm+=psiToXGammaAmp(1, 0, 0, -1, theData,PsiToEtaGamMag, PsiToEtaGamPhi );
  JmmGmp+=psiToXGammaAmp(-1, 0, 0, 1, theData, PsiToEtaGamMag, PsiToEtaGamPhi );
  JmmGmm+=psiToXGammaAmp(-1, 0, 0, -1, theData, PsiToEtaGamMag, PsiToEtaGamPhi  );
}


complex<double> JpsiGamEtaPiPiProdLh::A980DecFlatte(fitParams& theParamVal, const Vector4<double> &__p4){
  
  complex<double> result(0.,0.);

  const double massPi0 = 0.1349766;
  const double massEta = 0.547853; 
  std::pair <const double, const double> decPairPi0Eta=make_pair(massPi0, massEta);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  std::pair <const double, const double> decPairKK=make_pair(KplusMass,K0Mass);

  double a0_980Mass=theParamVal.Masses[paramEnumJpsiGamEtaPiPi::a0_980];
  double a0_980gPiEta=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gPiEta];
  double a0_980gKK=theParamVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gKK];

  result+=Flatte(__p4, decPairPi0Eta, decPairKK, a0_980Mass, a0_980gPiEta, a0_980gKK);

  return result;
}


void JpsiGamEtaPiPiProdLh::print(std::ostream& os) const{
  os << "JpsiGamEtaPiPiProdLh::print\n";
}


void JpsiGamEtaPiPiProdLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  JpsiGamEtaPiPiFitParams theFitParams;
  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;
  
  if(_etaToPiPiEtaHyp || _etaToa980PiHyp || _etaToa2_1320PiHyp){
    theAmpMap[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::PsiToEtaGamma);

    if(_etaToPiPiEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToPiPiEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToPiPiEta);
    }

    if (_etaToa980PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToA980Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToA980Pi);
    }
    
    if(_etaToa2_1320PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::EtaToA2_1320Pi);
    }
  }

  
  if(_f1ToPiPiEtaHyp || _f1Toa980PiHyp){
    theAmpMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::PsiToF1Gamma);
    
    if(_f1ToPiPiEtaHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1ToPiPiEta] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1ToPiPiEta);
    }

    if (_f1Toa980PiHyp){
      theAmpMap[paramEnumJpsiGamEtaPiPi::F1ToA980Pi] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::F1ToA980Pi);
    }
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
    fitVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gPiEta]= .30;
    fitErr.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gPiEta]=0.05;
    fitVal.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gKK]= .30;
    fitErr.gFactors[paramEnumJpsiGamEtaPiPi::a0_980gKK]=0.05;
  }

  if(_etaToa2_1320PiHyp){
    fitVal.Masses[paramEnumJpsiGamEtaPiPi::a2_1320]=1.32;
    fitErr.Masses[paramEnumJpsiGamEtaPiPi::a2_1320]=0.02;
    fitVal.Widths[paramEnumJpsiGamEtaPiPi::a2_1320]=0.1;
    fitErr.Widths[paramEnumJpsiGamEtaPiPi::a2_1320]=0.03;
  }
  //fill other params  
  if(_usePhasespace){
    fitVal.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace]=0.2;
    fitErr.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace]=0.4;
  }

}
  



bool 
JpsiGamEtaPiPiProdLh::initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   ){
  
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("etaToPiPiEtaHyp");
  
  if (iter !=hypMap.end()){
    _etaToPiPiEtaHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _etaToPiPiEtaHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis etaToPiPiEtaHyp not set!!!" <<endmsg;

  iter= hypMap.find("etaToa980PiHyp");
  if (iter !=hypMap.end()){
    _etaToa980PiHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _etaToa980PiHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis etaToa980PiHyp not set!!!" <<endmsg;    

  iter= hypMap.find("etaToa2_1320PiHyp");
  if (iter !=hypMap.end()){
    _etaToa2_1320PiHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _etaToa2_1320PiHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis etaToa2_1320PiHyp not set!!!" <<endmsg;    

  
  iter= hypMap.find("f1ToPiPiEtaHyp");
  if (iter !=hypMap.end()){
    _f1ToPiPiEtaHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f1ToPiPiEtaHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis f1ToPiPiEtaHyp not set!!!" <<endmsg;

  iter= hypMap.find("f1Toa980PiHyp");
  if (iter !=hypMap.end()){
    _f1Toa980PiHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f1Toa980PiHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis f1Toa980PiHyp not set!!!" <<endmsg;  

  iter= hypMap.find("usePhasespace");  
  if (iter !=hypMap.end()){
    _usePhasespace= iter->second;
    Info<< "Using phasespace for bg parameterization " << iter->first << "\t" << _usePhasespace <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "using phasespace not set!!!" <<endmsg;
  
  
  return true;
}








