#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKProdLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh(boost::shared_ptr<const EvtDataBaseList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsLh(theEvtList)
  ,_eta2225Hyp(true)
  ,_f02020Hyp(true)
  ,_f02020FlatteHyp(true)
  ,_f22300Hyp(true)
{
  
  
  initializeHypothesisMap( hypMap);
  


  
}

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsLh(theLhPtr->getEventList())
  ,_eta2225Hyp(true)
  ,_f02020Hyp(true)
  ,_f02020FlatteHyp(true)
  ,_f22300Hyp(true)
{
  
  initializeHypothesisMap( hypMap);
  
}

JpsiGamKsKlKKProdLh::~JpsiGamKsKlKKProdLh()
{;
}



double JpsiGamKsKlKKProdLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;

  Vector4<double> fvKsKlKpKm = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToEtacGamma];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToEtacGamma];
  
  double mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::etac];
  double width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::etac];
  
  complex<double> JmpGmp=etaGammaAmp(1, 0, 1, theData, PsiToEtacGamMag, PsiToEtacGamPhi, mass, width );
  complex<double> JmpGmm=etaGammaAmp(1, 0, -1, theData,  PsiToEtacGamMag, PsiToEtacGamPhi, mass, width );
  complex<double> JmmGmp=etaGammaAmp(-1, 0, 1, theData,  PsiToEtacGamMag, PsiToEtacGamPhi, mass, width );
  complex<double> JmmGmm=etaGammaAmp(-1, 0, -1, theData,  PsiToEtacGamMag, PsiToEtacGamPhi, mass, width );

  if (_eta2225Hyp){
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEta2225GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEta2225GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToEta2225Gamma];
    
    mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::eta2225];
    width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::eta2225];
    
    JmpGmp+=etaGammaAmp(1, 0, 1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi, mass, width );
    JmpGmm+=etaGammaAmp(1, 0, -1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi, mass, width );
    JmmGmp+=etaGammaAmp(-1, 0, 1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi, mass, width );
    JmmGmm+=etaGammaAmp(-1, 0, -1, theData, PsiToEta2225GamMag, PsiToEta2225GamPhi, mass, width );
  }

  if (_f02020Hyp || _f02020FlatteHyp ){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof02020GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof02020GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToF02020Gamma];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F02020ToPhiPhiMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F02020ToPhiPhiPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::F02020ToPhiPhi];
    
    
    complex<double> dynamicPart(0.,0.);
    Vector4<double> fv2Phi = theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
    mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::f02020];
    if( _f02020FlatteHyp  ){
      double gKK =  theParamVal.Widths[paramEnumJpsiGamKsKlKK::f02020gKK];
      double gPhiPhi =  theParamVal.Widths[paramEnumJpsiGamKsKlKK::f02020gPhiPhi];
      const double phiMass = 1.019455;
      const double kaonMass = 0.493677;
      std::pair <const double, const double> kkPair=make_pair(kaonMass, kaonMass);
      std::pair <const double, const double> phiphiPair=make_pair(phiMass, phiMass);
      dynamicPart = Flatte( fv2Phi , phiphiPair, kkPair,   mass, gPhiPhi, gKK  );
      
    } else{
      width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::f02020];
      dynamicPart = BreitWigner( fv2Phi, mass, width);
    }
    
    JmpGmp+=f0GammaAmp(1, 1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi ) * dynamicPart;
    JmpGmm+=f0GammaAmp(1, -1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi ) * dynamicPart;
    JmmGmp+=f0GammaAmp(-1, 1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi ) * dynamicPart;
    JmmGmm+=f0GammaAmp(-1, -1, theData,  PsiTof02020GamMag, PsiTof02020GamPhi,F02020ToPhiPhiMag,F02020ToPhiPhiPhi ) * dynamicPart;
  } 
  
  if (_f22300Hyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof2GamMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof2GamPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::PsiToF22300Gamma];
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F2ToPhiPhiMag=theParamVal.Mags[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > F2ToPhiPhiPhi=theParamVal.Phis[paramEnumJpsiGamKsKlKK::F22300ToPhiPhi];
    
    mass = theParamVal.Masses[paramEnumJpsiGamKsKlKK::f22300];
    width = theParamVal.Widths[paramEnumJpsiGamKsKlKK::f22300];
    
    JmpGmp+=f2GammaAmp(1, 1, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width );
    JmpGmp+=f2GammaAmp(1, -1, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width );
    JmpGmp+=f2GammaAmp(-1, 1, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width );
    JmpGmp+=f2GammaAmp(-1, -1, theData,  PsiTof2GamMag, PsiTof2GamPhi,F2ToPhiPhiMag,F2ToPhiPhiPhi,mass,width );
  } 
  
  
  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm); 
  return result;  
}

complex<double> JpsiGamKsKlKKProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}




complex<double> JpsiGamKsKlKKProdLh::etaGammaAmp(Spin Minit, Spin Metac, Spin Mgamma, EvtData* theData, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaMag, 
						 std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& PsiToEtaPhi,
						 double mass, double width){
  
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

   result*=BreitWigner( fv2Phi, mass, width)*etaToPhiPhiTo4KAmp(theData);

return result;
}







complex<double> JpsiGamKsKlKKProdLh::f0GammaAmp(Spin Minit, Spin Mgamma, EvtData* theData, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdMag, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0ProdPhi,
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecMag,  
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf0DecPhi ){
  
   complex<double> result(0.,0.);
   Spin f0Spin=0;
   //Vector4<double> fv2Phi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
   
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=ampf0ProdMag.begin(); itPsi!=ampf0ProdMag.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second;
     double thePsiPhi=ampf0ProdMag[PsiState];
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
   result*f0ToPhiPhiTo4KAmp(theData, ampf0DecMag,ampf0DecPhi );
   
   
return result;
}




complex<double> JpsiGamKsKlKKProdLh::etaToPhiPhiTo4KAmp(EvtData* theData){

  complex<double> result(0.,0.);

    for (Spin lambdaPhi=-1; lambdaPhi<1; lambdaPhi++){
       result+=Clebsch(1, lambdaPhi, 1, -lambdaPhi, 1, 0)*3.*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KsKl][1][lambdaPhi][0])*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KpKm][1][lambdaPhi][0]); //Clebsch(1, lambdaPhi, 1, -lambdaPhi, 1, 0): jj-coupling mother->phi phi; 3=sqrt(2L+1)*sqrt(2L+1) Cls=1


      }

    return result;
}

complex<double> JpsiGamKsKlKKProdLh::f0ToPhiPhiTo4KAmp(EvtData* theData, 
				    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf0DecMag , 
                 		    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf0DecPhi ){
  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itf0;
   for ( itf0=ampf0DecMag.begin(); itf0!=ampf0DecMag.end(); ++itf0){
     boost::shared_ptr<const JPCLS> f0ToPhiPhiState=itf0->first;
     double theMag=itf0->second;
     double thePhi=ampf0DecPhi[f0ToPhiPhiState];
     complex<double> expiphi(cos(thePhi), sin(thePhi));

     complex<double> tmp2PhiDecAmp(0.,0.);
     for (Spin lambdaPhi1=-1; lambdaPhi1<1; lambdaPhi1++){
       Spin lambdaPhi2=lambdaPhi1;
       tmp2PhiDecAmp+=Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, f0ToPhiPhiState->S, lambdaPhi1-lambdaPhi2)
	 *3.*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KsKl][1][lambdaPhi1][0])*conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KpKm][1][lambdaPhi2][0]); //3=sqrt(2L+1)*sqrt(2L+1) Cls=1    
     }
   
     result+= theMag*expiphi*sqrt(2*f0ToPhiPhiState->L+1)*tmp2PhiDecAmp;
   }
   return result;
}





complex<double> JpsiGamKsKlKKProdLh::f2GammaAmp(Spin Minit, Spin Mgamma, EvtData* theData, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdMag, 
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2ProdPhi,
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecMag,  
						std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ampf2DecPhi,
						double mass, double width ){
  
  complex<double> result(0.,0.);
  Vector4<double> fv2Phi= theData->FourVecs[enumJpsiGamKsKlKKData::V4_KsKlKpKm_HeliPsi];
  
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
  for ( itPsi=ampf2ProdMag.begin(); itPsi!=ampf2ProdMag.end(); ++itPsi){
    boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
    double thePsiMag=itPsi->second;
    double thePsiPhi=ampf2ProdMag[PsiState];
    complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
    
    for(Spin f2Lambda=-2;f2Lambda<=2; f2Lambda++){
      Spin lambda = f2Lambda-Mgamma;
      if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
      
      complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
	*Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
	*Clebsch(2, f2Lambda, 1, -Mgamma, PsiState->S, lambda  )
	*conj( theData->WignerDs[enumJpsiGamKsKlKKData::Df_Psi][PsiState->J][Minit][lambda]  );
      
      amp=amp*f2ToPhiPhiTo4KAmp(theData, f2Lambda,  ampf2DecMag,ampf2DecPhi );
      
      result+= amp;
    }
  }
  
  result*=BreitWigner( fv2Phi, mass, width);
  
  return result;
}


complex<double> JpsiGamKsKlKKProdLh::f2ToPhiPhiTo4KAmp( EvtData* theData, Spin f2Lambda, 
							std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecMag , 
							std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > &ampf2DecPhi  ){
  
  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itf2;
  for ( itf2=ampf2DecMag.begin(); itf2!=ampf2DecMag.end(); ++itf2){
    boost::shared_ptr<const JPCLS> f2State=itf2->first;
    double theMag=itf2->second;
    double thePhi=ampf2DecMag[f2State];
    complex<double> expiphi(cos(thePhi), sin(thePhi));
    
    for(Spin lambdaPhi1=-1; lambdaPhi1<=1; lambdaPhi1++){
      for(Spin lambdaPhi2=-1; lambdaPhi2<=1; lambdaPhi2++){
	Spin lambda = lambdaPhi1-lambdaPhi2;
	if( fabs(lambda)>f2State->J || fabs(lambda)>f2State->S) continue;
	
	complex<double> amp = theMag*expiphi*sqrt(2*f2State->L+1)
	  *Clebsch(f2State->L, 0, f2State->S, lambda, f2State->J, lambda)
	  *Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, f2State->S, lambda  )
	  *conj( theData->WignerDs[enumJpsiGamKsKlKKData::Df_Spin2][f2State->J][f2Lambda][lambda]  );
	
	amp = amp * phiphiTo4KAmp( theData, lambdaPhi1, lambdaPhi2 );
	result +=amp;
      }
    }
    
  }

  return result;
}


complex<double> JpsiGamKsKlKKProdLh::phiphiTo4KAmp( EvtData* theData, Spin lambdaPhi1, Spin lambdaPhi2  ){
  complex<double> result(0.,0.);
  
  result = 3. * conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KsKl][1][lambdaPhi1][0])
    * conj(theData->WignerDs[enumJpsiGamKsKlKKData::Df_KpKm][1][lambdaPhi2][0]);
  
  return result;
}





void JpsiGamKsKlKKProdLh::print(std::ostream& os) const{
  os << "JpsiGamKsKlKKProdLh::print\n";
}


void JpsiGamKsKlKKProdLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  JpsiGamKsKlKKFitParams theFitParams;
  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;
  theAmpMap[paramEnumJpsiGamKsKlKK::PsiToEtacGamma] = theFitParams.jpclsVec(paramEnumJpsiGamKsKlKK::PsiToEtacGamma);
  
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
  fitVal.Masses[paramEnumJpsiGamKsKlKK::etac]=2.95;
  fitErr.Masses[paramEnumJpsiGamKsKlKK::etac]=0.05;
  fitVal.Widths[paramEnumJpsiGamKsKlKK::etac]=.06;
  fitErr.Widths[paramEnumJpsiGamKsKlKK::etac]=0.03;
  
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
  
  return true;
}
