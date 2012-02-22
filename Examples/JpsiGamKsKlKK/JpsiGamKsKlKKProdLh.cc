#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKProdLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"
#include "Examples/JpsiGamKsKlKK/FitParamErrorMatrix.hh"
#include "Examples/JpsiGamKsKlKK/FitParamIndex.hh"

#include "ErrLogger/ErrLogger.hh"

#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh(boost::shared_ptr<const EvtDataBaseList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsLh(theEvtList)
  ,_etacHyp(false)
  ,_eta2225Hyp(false)
  ,_f02020Hyp(false)
  ,_f02020FlatteHyp(false)
  ,_f22300Hyp(false)
  ,_eta21870Hyp(false)
  ,_f1Hyp(false)
  ,_usePhasespace(false)
  ,_useCommonProductionPhase(true)
  ,_massIndependentFit(false)
  ,_phiMass( 1.019455)
  ,_kaonMass(0.493677)
{
  initializeHypothesisMap( hypMap);
 
}

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsLh(theLhPtr->getEventList())
  ,_etacHyp(false)
  ,_eta2225Hyp(false)
  ,_f02020Hyp(false)
  ,_f02020FlatteHyp(false)
  ,_f22300Hyp(false) 
  ,_eta21870Hyp(false)
  ,_f1Hyp(false)
  ,_usePhasespace(false)
  ,_useCommonProductionPhase(true)
  ,_massIndependentFit(false)
  ,_phiMass( 1.019455)
  ,_kaonMass(0.493677)
{
  
  initializeHypothesisMap( hypMap);
  
}

JpsiGamKsKlKKProdLh::~JpsiGamKsKlKKProdLh()
{;
}



double JpsiGamKsKlKKProdLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;
  
 
  complex<double> JmpGmp(0.,0.);
  complex<double> JmpGmm(0.,0.);
  complex<double> JmmGmp(0.,0.);
  complex<double> JmmGmm(0.,0.);
  
  if(_massIndependentFit){
    
    if(_etacHyp) calcEtacGammaAmp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::MassIndependent );
    if(_f02020Hyp) calcF02020Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent );
    if(_f22300Hyp) calcF22300Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent);
    if(_eta21870Hyp) calcE21870Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent);
    if(_f1Hyp) calcF1Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent);
    
    result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
    return result;
  }
  
  //
  //dynamic is parameterized by Breit-Wigner or Flatte
  //
  if(_etacHyp) calcEtacGammaAmp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::BreitWigner );
  if(_eta2225Hyp) calcEta2225Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::BreitWigner );
  if(_f02020Hyp) calcF02020Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner );
  if(_f02020FlatteHyp ) calcF02020Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::Flatte );
  if(_f22300Hyp) calcF22300Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner);
  if(_eta21870Hyp) calcE21870Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner);
  if(_f1Hyp) calcF1Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner);
  
  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  if(_usePhasespace){
    result = result + theParamVal.otherParams[paramEnumJpsiGamKsKlKK::phaseSpace];
  }
  
  return result;  
}

complex<double> JpsiGamKsKlKKProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}



complex<double> JpsiGamKsKlKKProdLh::etaGammaAmp(Spin Minit, Spin Metac, Spin Mgamma, EvtData* theData, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaMag, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaPhi ){
  
   complex<double> result(0.,0.);

   Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
   
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToEtaMag.begin(); itPsi!=PsiToEtaMag.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second;
     double thePsiPhi=PsiToEtaPhi[PsiState];
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = Metac-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(0, Metac, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->WignerDs[enumJpsiGamKsKlKKData::Df_Psi][PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }

   result*=etaToPhiPhiTo4KAmp(theData);

return result;
}







complex<double> JpsiGamKsKlKKProdLh::f0GammaAmp(Spin Minit, Spin Mgamma, EvtData* theData, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdMag, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdPhi,
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecMag,  
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecPhi,
						bool useCommonProductionPhase, dynamicModelParams& dynModPars){
  
   complex<double> result(0.,0.);
   Spin f0Spin=0;
   //Vector4<double> fv2Phi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
   
   double theCommonPhasePhi=999.;
   bool setCommonPhase = false;
   
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=ampf0ProdMag.begin(); itPsi!=ampf0ProdMag.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second;
     double thePsiPhi=ampf0ProdPhi[PsiState];
     
     if(useCommonProductionPhase){
       if(!setCommonPhase){
	 theCommonPhasePhi=thePsiPhi;
	 setCommonPhase=true;
       }else{
	 thePsiPhi = theCommonPhasePhi;
       }
     }
     
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = f0Spin-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(0, f0Spin, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->WignerDs[enumJpsiGamKsKlKKData::Df_Psi][PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }

   //result*=BreitWigner( fv2Phi, mass, width)*f0ToPhiPhiTo4KAmp(theData, ampf0DecMag,ampf0DecPhi );
   result*=f0ToPhiPhiTo4KAmp(theData, ampf0DecMag,ampf0DecPhi, dynModPars );
   
   
return result;
}




complex<double> JpsiGamKsKlKKProdLh::etaToPhiPhiTo4KAmp(EvtData* theData){

  complex<double> result(0.,0.);

    for (Spin lambdaPhi=-1; lambdaPhi<=1; lambdaPhi++){
       result+=Clebsch(1, lambdaPhi, 1, -lambdaPhi, 1, 0)*3.*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KsKl][1][lambdaPhi][0])*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KpKm][1][lambdaPhi][0]); //Clebsch(1, lambdaPhi, 1, -lambdaPhi, 1, 0): jj-coupling mother->phi phi; 3=sqrt(2L+1)*sqrt(2L+1) Cls=1


      }

    return result;
}

complex<double> JpsiGamKsKlKKProdLh::f0ToPhiPhiTo4KAmp(EvtData* theData, 
						       std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf0DecMag , 
						       std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf0DecPhi,
						       dynamicModelParams& dynModPars){
  complex<double> result(0.,0.);
  complex<double> dynamicPart(0.,0.);
  Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  if(dynModPars.dynamicModel == dynamicModelParams::Flatte || dynModPars.dynamicModel==dynamicModelParams::BreitWigner){
    dynamicPart = dynModPars.value;
  }else if( dynModPars.dynamicModel==dynamicModelParams::MassIndependent  ){
    dynamicPart = complex<double>(1,0); 
  }
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itf0;
   for ( itf0=ampf0DecMag.begin(); itf0!=ampf0DecMag.end(); ++itf0){
     boost::shared_ptr<const JPCLS> f0ToPhiPhiState=itf0->first;
     double theMag=itf0->second;
     double thePhi=ampf0DecPhi[f0ToPhiPhiState];
     complex<double> expiphi(cos(thePhi), sin(thePhi));
     
     complex<double> tmp2PhiDecAmp(0.,0.);
     for (Spin lambdaPhi1=-1; lambdaPhi1<=1; lambdaPhi1++){
       Spin lambdaPhi2=lambdaPhi1;
       tmp2PhiDecAmp+=Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, f0ToPhiPhiState->S, lambdaPhi1-lambdaPhi2)
	 *3.*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KsKl][1][lambdaPhi1][0])*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KpKm][1][lambdaPhi2][0]); //3=sqrt(2L+1)*sqrt(2L+1) Cls=1    
     }
     
     if(dynModPars.dynamicModel==dynamicModelParams::BreitWignerBlattW){
       dynamicPart=BreitWignerBlattW(fv2Phi, _phiMass, _phiMass,  dynModPars.mass, dynModPars.width, f0ToPhiPhiState->L );
     }
     
     result+= theMag*expiphi*sqrt(2*f0ToPhiPhiState->L+1)*tmp2PhiDecAmp*dynamicPart;
   }
   return result;
}




complex<double> JpsiGamKsKlKKProdLh::f2GammaAmp(Spin Minit, Spin Mgamma, Spin fJSpin, EvtData* theData, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdMag, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdPhi,
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecMag,  
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecPhi,
						double mass, double width, bool useCommonProductionPhase, dynamicModelParams& dynModPars){
  
  complex<double> result(0.,0.);
  Vector4<double> fv2Phi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];

  std::map<Spin, complex<double> > decAmp;  
  for (Spin lambdaf2=-fJSpin; lambdaf2<=fJSpin; ++lambdaf2){
    decAmp[lambdaf2]=f2ToPhiPhiTo4KAmp(theData, lambdaf2,  ampf2DecMag,ampf2DecPhi, dynModPars );
  }
  
  double theCommonPhasePhi=999.;
  bool setCommonPhase=false;
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
  for ( itPsi=ampf2ProdMag.begin(); itPsi!=ampf2ProdMag.end(); ++itPsi){
    boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
    double thePsiMag=itPsi->second;
    double thePsiPhi=ampf2ProdPhi[PsiState];
    
    if(useCommonProductionPhase){
      if(!setCommonPhase){
	theCommonPhasePhi=thePsiPhi;
	setCommonPhase=true;
      }else{
	thePsiPhi = theCommonPhasePhi;
      }
    }
    
    complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
    
    for(Spin f2Lambda=-fJSpin;f2Lambda<=fJSpin; f2Lambda++){
      Spin lambda = f2Lambda-Mgamma;
      if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
      
      complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
	*Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
	*Clebsch(fJSpin, f2Lambda, 1, -Mgamma, PsiState->S, lambda  )
	*conj( theData->WignerDs[enumJpsiGamKsKlKKData::Df_Psi][PsiState->J][Minit][lambda]  );
      
      amp=amp*decAmp[f2Lambda];
      
      result+= amp;
    }
  }
  
  
  return result;
}


complex<double> JpsiGamKsKlKKProdLh::f2ToPhiPhiTo4KAmp( EvtData* theData, Spin f2Lambda, 
							std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecMag , 
							std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecPhi,
							dynamicModelParams& dynModPars  ){
  
  complex<double> result(0.,0.);
  
  complex<double> dynamicPart(0.,0.);
  Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  if(dynModPars.dynamicModel == dynamicModelParams::Flatte || dynModPars.dynamicModel==dynamicModelParams::BreitWigner){
    dynamicPart = dynModPars.value;
  }else if( dynModPars.dynamicModel==dynamicModelParams::MassIndependent  ){
    dynamicPart = complex<double>(1,0); 
  }
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itf2;
  for ( itf2=ampf2DecMag.begin(); itf2!=ampf2DecMag.end(); ++itf2){
    boost::shared_ptr<const JPCLS> f2State=itf2->first;
    double theMag=itf2->second;
    double thePhi=ampf2DecPhi[f2State];
    complex<double> expiphi(cos(thePhi), sin(thePhi));
    
    for(Spin lambdaPhi1=-1; lambdaPhi1<=1; lambdaPhi1++){
      for(Spin lambdaPhi2=-1; lambdaPhi2<=1; lambdaPhi2++){
	Spin lambda = lambdaPhi1-lambdaPhi2;
	if( fabs(lambda)>f2State->J || fabs(lambda)>f2State->S) continue;
	
	complex<double> amp = theMag*expiphi*sqrt(2*f2State->L+1)
	  *Clebsch(f2State->L, 0, f2State->S, lambda, f2State->J, lambda)
	  *Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, f2State->S, lambda  )
	  *conj( theData->WignerDs[enumJpsiGamKsKlKKData::Df_Spin2][f2State->J][f2Lambda][lambda]  );
	
	if(dynModPars.dynamicModel==dynamicModelParams::BreitWignerBlattW){
	  dynamicPart=BreitWignerBlattW(fv2Phi, _phiMass, _phiMass,  dynModPars.mass, dynModPars.width, f2State->L );
	}
	
 	amp = amp * phiphiTo4KAmp( theData, lambdaPhi1, lambdaPhi2 ) * dynamicPart;

	result +=amp;
      }
    }
    
  }

  return result;
}


complex<double> JpsiGamKsKlKKProdLh::phiphiTo4KAmp( EvtData* theData, Spin lambdaPhi1, Spin lambdaPhi2  ){
  complex<double> result(0.,0.);
  
  result = 3. * conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KsKl][1][lambdaPhi1][0])
    * 3.* conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KpKm][1][lambdaPhi2][0]);
  
  return result;
}




bool  JpsiGamKsKlKKProdLh::initializeFlatteModel( dynamicModelParams &theDynModel, const Vector4<double> &fv2Phi, double mass, double gPhiPhi, double gKK   ){
  
  theDynModel.dynamicModel = dynamicModelParams::Flatte;
  theDynModel.gFactor1 = gPhiPhi;
  theDynModel.gFactor2 = gKK;
  pair<const double, const double> mp1 = make_pair(_phiMass, _phiMass);
  pair<const double, const double> mp2 = make_pair(_kaonMass, _kaonMass);
  theDynModel.massPair1 = mp1; 
  theDynModel.massPair2 = mp2; 
  theDynModel.value = Flatte( fv2Phi, mp1, mp2, theDynModel.mass, theDynModel.gFactor1, theDynModel.gFactor2  );
  
  return true;
}


bool  JpsiGamKsKlKKProdLh::initializeBreitWignerModel(dynamicModelParams &theDynModel,const Vector4<double> &fv2Phi, double mass, double width, dynamicModelParams::enumDynamicModel theModel ){
  theDynModel.mass = mass;
  theDynModel.width=width;
  theDynModel.dynamicModel = theModel;
  if(theModel == dynamicModelParams::BreitWigner){
    theDynModel.value = BreitWigner(fv2Phi, mass, width );
  }
  return true;
}

    












void JpsiGamKsKlKKProdLh::print(std::ostream& os) const{
  os << "JpsiGamKsKlKKProdLh::print\n";
}


void JpsiGamKsKlKKProdLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  JpsiGamKsKlKKFitParams theFitParams;
  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;
  
  if(_etacHyp){
    theAmpMap[paramEnumJpsiGamKsKlKK::PsiToEtacGamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToEtacGamma);
  }
  
  if(_eta2225Hyp){
    theAmpMap[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma);
  }
  
  if(_f02020Hyp || _f02020FlatteHyp ){
    theAmpMap[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToF02020Gamma);
    theAmpMap[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::F02020ToPhiPhi);
  }

  if(_f22300Hyp){
    theAmpMap[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToF22300Gamma);
    theAmpMap[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::F22300ToPhiPhi);
  }
  
  if(_eta21870Hyp){
    theAmpMap[paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma);
    theAmpMap[paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi);
  }
  
  if(_f1Hyp){
    theAmpMap[paramEnumJpsiGamKsKlKK::PsiToF1Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToF1Gamma);
    theAmpMap[paramEnumJpsiGamKsKlKK::F1ToPhiPhi] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::F1ToPhiPhi);
  }
  
  
  
  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > >::iterator itAmpMap;
  for (itAmpMap=theAmpMap.begin(); itAmpMap!=theAmpMap.end(); ++itAmpMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > valMagMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > errMagMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > valPhiMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > errPhiMap;
    
    std::vector< boost::shared_ptr<const JPCLS> >::iterator itAmp;
    for (itAmp=itAmpMap->second.begin(); itAmp!=itAmpMap->second.end(); ++itAmp){
      valMagMap[(*itAmp)]=0.6;
      errMagMap[(*itAmp)]=0.3; 
      valPhiMap[(*itAmp)]=0.;
      errPhiMap[(*itAmp)]=0.2;      
    }
    
    fitVal.Mags[itAmpMap->first]=valMagMap;
    fitVal.Phis[itAmpMap->first]=valPhiMap;  
    fitErr.Mags[itAmpMap->first]=errMagMap;
    fitErr.Phis[itAmpMap->first]=errPhiMap;  
  }

  
  //fill masses and wisths
  if(_etacHyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::etac]=2.98;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::etac]=0.03;
    fitVal.Widths[paramEnumJpsiGamKsKlKK::etac]=.06;
    fitErr.Widths[paramEnumJpsiGamKsKlKK::etac]=0.03;
  }
  
  if(_eta2225Hyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::eta2225]=2.225;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::eta2225]=0.3;
    fitVal.Widths[paramEnumJpsiGamKsKlKK::eta2225]=.3;
    fitErr.Widths[paramEnumJpsiGamKsKlKK::eta2225]=0.1;
  }
  
  if(_f02020Hyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::f02020]=2.020;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::f02020]=0.3;
    fitVal.Widths[paramEnumJpsiGamKsKlKK::f02020]=.30;
    fitErr.Widths[paramEnumJpsiGamKsKlKK::f02020]=0.1;
  }

  if(_f02020FlatteHyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::f02020]=2.020;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::f02020]=0.3;
    fitVal.gFactors[paramEnumJpsiGamKsKlKK::f02020gKK]=.30;
    fitErr.gFactors[paramEnumJpsiGamKsKlKK::f02020gKK]=0.1;
    fitVal.gFactors[paramEnumJpsiGamKsKlKK::f02020gPhiPhi]=.7;
    fitErr.gFactors[paramEnumJpsiGamKsKlKK::f02020gPhiPhi]=0.3;
  }
  
  
  if(_f22300Hyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::f22300]=2.186;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::f22300]=0.3;
    fitVal.Widths[paramEnumJpsiGamKsKlKK::f22300]=.55;
    fitErr.Widths[paramEnumJpsiGamKsKlKK::f22300]=0.4;
  }
  
  if(_eta21870Hyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::eta21870]=2.1;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::eta21870]=0.3;
    fitVal.Widths[paramEnumJpsiGamKsKlKK::eta21870]=.55;
    fitErr.Widths[paramEnumJpsiGamKsKlKK::eta21870]=0.4;
  }
  
  if(_f1Hyp){
    fitVal.Masses[paramEnumJpsiGamKsKlKK::f1]=2.1;
    fitErr.Masses[paramEnumJpsiGamKsKlKK::f1]=0.3;
    fitVal.Widths[paramEnumJpsiGamKsKlKK::f1]=.55;
    fitErr.Widths[paramEnumJpsiGamKsKlKK::f1]=0.4;
  }
  
  if(_usePhasespace){
    fitVal.otherParams[paramEnumJpsiGamKsKlKK::phaseSpace]=0.2;
    fitErr.otherParams[paramEnumJpsiGamKsKlKK::phaseSpace]=0.4;
  }

}
  



bool 
JpsiGamKsKlKKProdLh::initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   ){
  
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("eta2225Hyp");
  
  if (iter !=hypMap.end()){
    _eta2225Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _eta2225Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis eta2225Hyp not set!!!" <<endmsg;
  
  
  iter= hypMap.find("etacHyp");
  if (iter !=hypMap.end()){
    _etacHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _etacHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis etacHyp not set!!!" <<endmsg;
  

  iter= hypMap.find("f02020Hyp");  
  if (iter !=hypMap.end()){
    _f02020Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f02020Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis f02020Hyp not set!!!" <<endmsg;

  iter= hypMap.find("f02020FlatteHyp");  
  if (iter !=hypMap.end()){
    _f02020FlatteHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f02020FlatteHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis f02020FlatteHyp not set!!!" <<endmsg;
  
  
  iter= hypMap.find("f22300Hyp");  
  if (iter !=hypMap.end()){
    _f22300Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f22300Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis f22300Hyp not set!!!" <<endmsg;
  
  iter= hypMap.find("eta21870Hyp");  
  if (iter !=hypMap.end()){
    _eta21870Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _eta21870Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis eta21870Hyp not set!!!" <<endmsg;

  iter= hypMap.find("f1Hyp");  
  if (iter !=hypMap.end()){
    _f1Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f1Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis f1Hyp not set!!!" <<endmsg;
  

  iter= hypMap.find("usePhasespace");  
  if (iter !=hypMap.end()){
    _usePhasespace= iter->second;
    Info<< "Using phasespace for bg parameterization " << iter->first << "\t" << _usePhasespace <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "using phasespace not set!!!" <<endmsg;
  
  
  return true;
}


bool 
JpsiGamKsKlKKProdLh::calcComponentIntensity(  EvtData* theData, fitParams& theParamVal, FitParamErrorMatrix& theErrMatrix, std::string component, std::pair<double, double> &intensity    ){
  
  complex<double> JmpGmp(0.0,0.0);
  complex<double> JmpGmm(0.0,0.0);
  complex<double> JmmGmp(0.0,0.0);
  complex<double> JmmGmm(0.0,0.0);
  
  intensity = make_pair(0.,0.);
  
  
  std::map<const std::string, bool>::const_iterator iter= _hypMap.find(component);
  if (iter ==_hypMap.end()){
    static int alerts=0;
    if(alerts<5) Alert << "Component " << component << " was not included in fit hypothesis";
    alerts++;
    return true;
  }
  if (iter !=_hypMap.end() && !iter->second ){
    static int alerts=0;
    if(alerts<5) Alert << "Component " << component << " was disabled in fit hypothesis" ;
    alerts++;
    return true;
  }
  
  
  if(_massIndependentFit){
    
    if(component == "etacHyp" ){
      calcEtacGammaAmp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::MassIndependent  );
      intensity.first =norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
      
      boost::function<void(EvtData* , fitParams&, complex<double>&, complex<double>&, complex<double>&, complex<double>&, dynamicModelParams::enumDynamicModel) > ampFunction = boost::bind(&JpsiGamKsKlKKProdLh::calcEtacGammaAmp, this, _1, _2, _3, _4, _5, _6,_7  );
      double intensityError=0.0;
      std::vector< int > theAmps;
      theAmps.push_back( paramEnumJpsiGamKsKlKK::PsiToEtacGamma );
      calcAmpError(  theData, theParamVal, theErrMatrix, dynamicModelParams::MassIndependent, ampFunction, theAmps, intensityError  );
      intensity.second=intensityError;
    }
    
    
    else if( component == "f02020Hyp" ){
      calcF02020Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent );
      intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
      boost::function<void(EvtData* , fitParams&, complex<double>&, complex<double>&, complex<double>&, complex<double>&, dynamicModelParams::enumDynamicModel) > ampFunction = boost::bind(&JpsiGamKsKlKKProdLh::calcF02020Amp, this, _1, _2, _3, _4, _5, _6,_7  );
      double intensityError=0.0;
      std::vector< int > theAmps;
      theAmps.push_back( paramEnumJpsiGamKsKlKK::PsiToF02020Gamma );
      theAmps.push_back( paramEnumJpsiGamKsKlKK::F02020ToPhiPhi );
      calcAmpError(  theData, theParamVal, theErrMatrix, dynamicModelParams::MassIndependent, ampFunction, theAmps, intensityError  );
      intensity.second=intensityError;
    }
    
    else if(component == "f22300Hyp"){
      calcF22300Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent);
      intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
      boost::function<void(EvtData* , fitParams&, complex<double>&, complex<double>&, complex<double>&, complex<double>&, dynamicModelParams::enumDynamicModel) > ampFunction = boost::bind(&JpsiGamKsKlKKProdLh::calcF22300Amp, this, _1, _2, _3, _4, _5, _6,_7  );
      double intensityError=0.0;
      std::vector< int > theAmps;
      theAmps.push_back( paramEnumJpsiGamKsKlKK::PsiToF22300Gamma );
      theAmps.push_back( paramEnumJpsiGamKsKlKK::F22300ToPhiPhi );
      calcAmpError(  theData, theParamVal, theErrMatrix, dynamicModelParams::MassIndependent, ampFunction, theAmps, intensityError  );
      intensity.second=intensityError;
    }
    
    else if(component == "eta21870Hyp"){
      calcE21870Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent);
      intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
      boost::function<void(EvtData* , fitParams&, complex<double>&, complex<double>&, complex<double>&, complex<double>&, dynamicModelParams::enumDynamicModel) > ampFunction = boost::bind(&JpsiGamKsKlKKProdLh::calcE21870Amp, this, _1, _2, _3, _4, _5, _6,_7  );
      double intensityError=0.0;
      std::vector< int > theAmps;
      theAmps.push_back( paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma );
      theAmps.push_back( paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi );
      calcAmpError(  theData, theParamVal, theErrMatrix, dynamicModelParams::MassIndependent, ampFunction, theAmps, intensityError  );
      intensity.second=intensityError;
      
    }
    
    else if(component == "f1Hyp"){
      calcF1Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::MassIndependent);
      intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
    
      boost::function<void(EvtData* , fitParams&, complex<double>&, complex<double>&, complex<double>&, complex<double>&, dynamicModelParams::enumDynamicModel) > ampFunction = boost::bind(&JpsiGamKsKlKKProdLh::calcF1Amp, this, _1, _2, _3, _4, _5, _6,_7  );
      double intensityError=0.0;
      std::vector< int > theAmps;
      theAmps.push_back( paramEnumJpsiGamKsKlKK::PsiToF1Gamma );
      theAmps.push_back( paramEnumJpsiGamKsKlKK::F1ToPhiPhi );
      calcAmpError(  theData, theParamVal, theErrMatrix, dynamicModelParams::MassIndependent, ampFunction, theAmps, intensityError  );
      intensity.second=intensityError;


    }
    return true;
  }
  
  
  
  //
  //mass-dependent parameterization
  //
  
  if(component == "etacHyp" ){
    calcEtacGammaAmp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner    );
    intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  }
  
  else if(component == "eta2225Hyp"){
    calcEta2225Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::BreitWigner );
  intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  }


  else if( component == "f02020Hyp" ){
    calcF02020Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner );
  intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  }
  
  else if(component == "f02020FlatteHyp"){
    calcF02020Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::Flatte );
  intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  }
  
  else if(component == "f22300Hyp"){
    calcF22300Amp( theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner);
    intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  }
  
  else if(component == "eta21870Hyp"){
    calcE21870Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner);
    intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);  
}
  
  else if(component == "f1Hyp"){
    calcF1Amp(theData, theParamVal, JmpGmp,JmpGmm, JmmGmp,JmmGmm, dynamicModelParams::BreitWigner);
    intensity.first=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
  }
  
  return  true;
}


 void 
 JpsiGamKsKlKKProdLh::calcEtacGammaAmp( EvtData* theData, fitParams& theParamVal, 
					complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
					dynamicModelParams::enumDynamicModel theModel){
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToEtacGamma];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToEtacGamma];
   

   Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
   double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::etac];
   double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::etac];
   complex<double> theDynamicPart(1.,0.);
   if(theModel == dynamicModelParams::BreitWigner) theDynamicPart = BreitWigner(fv2Phi, mass, width);
   //complex<double> theDynamicPart = BreitWignerBlattW(fv2Phi, _phiMass, _phiMass, mass, width, 1);
   
   JmpGmp+=etaGammaAmp(1, 0, 1, theData, PsiToEtacGamMag, PsiToEtacGamPhi )*theDynamicPart;
   JmpGmm+=etaGammaAmp(1, 0, -1, theData,  PsiToEtacGamMag, PsiToEtacGamPhi  )*theDynamicPart;
   JmmGmp+=etaGammaAmp(-1, 0, 1, theData,  PsiToEtacGamMag, PsiToEtacGamPhi )*theDynamicPart;
   JmmGmm+=etaGammaAmp(-1, 0, -1, theData,  PsiToEtacGamMag, PsiToEtacGamPhi  )*theDynamicPart;
 }



void 
JpsiGamKsKlKKProdLh::calcEta2225Amp(EvtData* theData, fitParams& theParamVal, 
				    complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm,
				    dynamicModelParams::enumDynamicModel theModel){
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEta2225GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEta2225GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma];
  
  Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::eta2225];
  double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::eta2225];
  complex<double> theDynamicPart(1.,0.);
  if(theModel == dynamicModelParams::BreitWigner) theDynamicPart = BreitWigner(fv2Phi, mass, width);
  //complex<double> theDynamicPart = BreitWignerBlattW(fv2Phi, _phiMass, _phiMass, mass, width, 1);
  
  JmpGmp+=etaGammaAmp(1, 0, 1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi )*theDynamicPart;
  JmpGmm+=etaGammaAmp(1, 0, -1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi )*theDynamicPart;
  JmmGmp+=etaGammaAmp(-1, 0, 1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi )*theDynamicPart;
  JmmGmm+=etaGammaAmp(-1, 0, -1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi )*theDynamicPart;
  
}



void 
JpsiGamKsKlKKProdLh::calcF02020Amp(EvtData* theData, fitParams& theParamVal, 
				   complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
				   dynamicModelParams::enumDynamicModel theModel){
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof02020GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof02020GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma];
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F02020ToPhiPhiMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F02020ToPhiPhiPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi];
  
  Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  dynamicModelParams theDynModel;
  double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::f02020];
  
  if( theModel == dynamicModelParams::Flatte  ){
    double gKK =  theParamVal.gFactors[paramEnumJpsiGamKsKlKK::f02020gKK];
    double gPhiPhi =  theParamVal.gFactors[paramEnumJpsiGamKsKlKK::f02020gPhiPhi];
    initializeFlatteModel( theDynModel, fv2Phi, mass, gPhiPhi, gKK   );
    
  }else if( theModel == dynamicModelParams::BreitWigner ){
    double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::f02020];
    initializeBreitWignerModel( theDynModel, fv2Phi,  mass, width, dynamicModelParams::BreitWigner  );
  }else if( theModel == dynamicModelParams::MassIndependent ){
    theDynModel.value=complex<double>(1.,0.);
    theDynModel.dynamicModel = dynamicModelParams::MassIndependent;
  }

  
  
  JmpGmp+=f0GammaAmp(1, 1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi,_useCommonProductionPhase, theDynModel );
  JmpGmm+=f0GammaAmp(1, -1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi,_useCommonProductionPhase, theDynModel );
  JmmGmp+=f0GammaAmp(-1, 1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi,_useCommonProductionPhase, theDynModel );
  JmmGmm+=f0GammaAmp(-1, -1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi,_useCommonProductionPhase, theDynModel );
}



void JpsiGamKsKlKKProdLh::calcF22300Amp(EvtData* theData, fitParams& theParamVal, 
					complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
					dynamicModelParams::enumDynamicModel theModel){
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof2GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof2GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma];
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F2ToPhiPhiMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F2ToPhiPhiPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi];
  
  Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::f22300];
  double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::f22300];
  
  Spin f2Spin=2;
  dynamicModelParams theDynModel;
  if(theModel==dynamicModelParams::BreitWigner) {
    initializeBreitWignerModel( theDynModel, fv2Phi,  mass, width, dynamicModelParams::BreitWigner  ); 
  }
  else if(theModel == dynamicModelParams::MassIndependent){
    theDynModel.value=complex<double>(1.,0);
    theDynModel.dynamicModel=dynamicModelParams::MassIndependent;
  }
  
  
  
  JmpGmp+=f2GammaAmp(1, 1, f2Spin, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
  JmpGmm+=f2GammaAmp(1, -1, f2Spin, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
  JmmGmp+=f2GammaAmp(-1, 1, f2Spin, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
  JmmGmm+=f2GammaAmp(-1, -1, f2Spin, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
} 

 



void 
JpsiGamKsKlKKProdLh::calcE21870Amp(EvtData* theData, fitParams& theParamVal, 
				   complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
				   dynamicModelParams::enumDynamicModel theModel){
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEta2GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEta2GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToEta21870Gamma];
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > Eta2ToPhiPhiMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > Eta2ToPhiPhiPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::Eta21870ToPhiPhi];
  
  Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::eta21870];
  double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::eta21870];
  
  Spin eta2Spin=2;
  dynamicModelParams theDynModel;
  if(theModel==dynamicModelParams::BreitWigner){ initializeBreitWignerModel( theDynModel, fv2Phi,  mass, width, dynamicModelParams::BreitWigner  );
  }
  else if(theModel == dynamicModelParams::MassIndependent){
    theDynModel.value=complex<double>(1.,0);
    theDynModel.dynamicModel=dynamicModelParams::MassIndependent;
  }
  
  JmpGmp+=f2GammaAmp(1, 1, eta2Spin, theData,  PsiToEta2GamMag, PsiToEta2GamPhi,Eta2ToPhiPhiMag,Eta2ToPhiPhiPhi,mass,width, _useCommonProductionPhase,theDynModel );
  JmpGmm+=f2GammaAmp(1, -1, eta2Spin,theData,  PsiToEta2GamMag, PsiToEta2GamPhi,Eta2ToPhiPhiMag,Eta2ToPhiPhiPhi,mass,width, _useCommonProductionPhase,theDynModel );
  JmmGmp+=f2GammaAmp(-1, 1, eta2Spin,theData,  PsiToEta2GamMag, PsiToEta2GamPhi,Eta2ToPhiPhiMag,Eta2ToPhiPhiPhi,mass,width, _useCommonProductionPhase,theDynModel );
  JmmGmm+=f2GammaAmp(-1, -1, eta2Spin, theData,  PsiToEta2GamMag, PsiToEta2GamPhi,Eta2ToPhiPhiMag,Eta2ToPhiPhiPhi,mass,width , _useCommonProductionPhase,theDynModel );
}     




 void 
 JpsiGamKsKlKKProdLh::calcF1Amp(EvtData* theData, fitParams& theParamVal, 
				complex<double> &JmpGmp, complex<double> &JmpGmm, complex<double> &JmmGmp, complex<double> &JmmGmm, 
				dynamicModelParams::enumDynamicModel theModel){
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof1GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToF1Gamma];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof1GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToF1Gamma];
   
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F1ToPhiPhiMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::F1ToPhiPhi];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F1ToPhiPhiPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::F1ToPhiPhi];
   
   Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
   double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::f1];
   double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::f1];
   
   Spin f1Spin=1;
   dynamicModelParams theDynModel;
   if(theModel == dynamicModelParams::BreitWigner ){
     initializeBreitWignerModel( theDynModel, fv2Phi,  mass, width, dynamicModelParams::BreitWigner  );
   }
   else if(theModel == dynamicModelParams::MassIndependent){
     theDynModel.value=complex<double>(1.,0);
     theDynModel.dynamicModel=dynamicModelParams::MassIndependent;
   }
   
   JmpGmp+=f2GammaAmp(1, 1, f1Spin, theData,  PsiTof1GamMag, PsiTof1GamPhi,F1ToPhiPhiMag,F1ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
   JmpGmm+=f2GammaAmp(1, -1, f1Spin, theData,  PsiTof1GamMag, PsiTof1GamPhi,F1ToPhiPhiMag,F1ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
   JmmGmp+=f2GammaAmp(-1, 1, f1Spin, theData,  PsiTof1GamMag, PsiTof1GamPhi,F1ToPhiPhiMag,F1ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
   JmmGmm+=f2GammaAmp(-1, -1, f1Spin, theData,  PsiTof1GamMag, PsiTof1GamPhi,F1ToPhiPhiMag,F1ToPhiPhiPhi,mass,width, _useCommonProductionPhase, theDynModel );
 }



 void
 JpsiGamKsKlKKProdLh::dumpComponentIntensity( std::ostream &os, fitParams& theParams, FitParamErrorMatrix &theErrMatrix  ){
   
   const std::vector<EvtData*> mcList=getEventList()->getMcVecs();
   const std::vector<EvtData*> dataList=getEventList()->getDataVecs();
   
   std::map<const std::string, bool>::const_iterator hypo= _hypMap.begin();
   while(hypo !=_hypMap.end()){
     if( hypo->second ){
       
       std::vector<EvtData*>::const_iterator it=mcList.begin();
       double integral=0.0;
       double integralError=0.0;
       while(it!=mcList.end()){
	 
	 std::pair<double, double> intensityEvent=make_pair(0.,0.);
	 calcComponentIntensity( *it, theParams, theErrMatrix, hypo->first, intensityEvent );
	 integral+=intensityEvent.first;
	 integralError+=intensityEvent.second;
	 it++;
       }
       Info << "Unscaled intensity for component " << hypo->first << ": " << integral << " +/- " << integralError << endmsg;
       double scale = (1.*dataList.size()) / (1.*mcList.size()) ;
       Info << "Scale factor: " << scale << endmsg;
       integral*=scale;
       integralError*=scale;
       Info << "Events for component " << hypo->first << ": " << integral << " +/- " << integralError << endmsg;
       
       os  << hypo->first << "\t" << integral << endl;
     }
     hypo++;
   }
 }





void JpsiGamKsKlKKProdLh::calcAmpError( EvtData* theData, fitParams& theParamVal, FitParamErrorMatrix& theErrMatrix,  
					dynamicModelParams::enumDynamicModel theModel, 
					boost::function<void(EvtData* , fitParams&, complex<double>&, 
							     complex<double>&, complex<double>&, complex<double>&, 
							     dynamicModelParams::enumDynamicModel)> calcAmp, 
					std::vector< int > theAmpsEnum,
					double& theAmpError)
{
  
  double epsilon=0.5e-5;
  double error(0.);
  std::vector<double> derivatives;
  std::vector<int> parameterIndices;
  
  std::vector<int>::iterator ampIter;
  for( ampIter=theAmpsEnum.begin(); ampIter!=theAmpsEnum.end(); ampIter++){
    
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > theMag = theParamVal.Mags[*ampIter];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator it;
    for ( it=theMag.begin(); it!=theMag.end(); ++it){
      
      fitParams newParamsPlus = theParamVal;
      fitParams newParamsMinus = theParamVal;
      boost::shared_ptr<const JPCLS> theState=it->first;
      newParamsPlus.Mags[*ampIter][theState]=theParamVal.Mags[*ampIter][theState]+epsilon;
      newParamsMinus.Mags[*ampIter][theState]=theParamVal.Mags[*ampIter][theState]-epsilon;
      
      complex<double> JmpGmp(0.0,0.0);
      complex<double> JmpGmm(0.0,0.0);
      complex<double> JmmGmp(0.0,0.0);
      complex<double> JmmGmm(0.0,0.0);
    
      calcAmp( theData, newParamsPlus, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::MassIndependent  );
      double resultPlus=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
      
      JmpGmp = complex<double>(0.0,0.0);
      JmpGmm =  complex<double>(0.0,0.0);
      JmmGmp =  complex<double>(0.0,0.0);
      JmmGmm =  complex<double>(0.0,0.0);
      
      calcAmp( theData, newParamsMinus, JmpGmp,JmpGmm, JmmGmp,JmmGmm,dynamicModelParams::MassIndependent  );
      double resultMinus=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm);
      
      derivatives.push_back((resultPlus-resultMinus)/(2*epsilon));
      static FitParamIndex theParamIndex(theParamVal);
      parameterIndices.push_back( theParamIndex.Mag(*ampIter, theState ) );
    }
    
    int elements = derivatives.size();
    for(int i=0; i<elements; i++){
      for(int j=i;j<elements;j++){
	if(i==j){
	  error+=derivatives[i]*derivatives[j]*theErrMatrix(parameterIndices[i],parameterIndices[j] );
	}else{
	  error+=derivatives[i]*derivatives[j]*theErrMatrix(parameterIndices[i],parameterIndices[j] );
	}
      }
    }
  }    
  theAmpError=sqrt(error);
  
}



