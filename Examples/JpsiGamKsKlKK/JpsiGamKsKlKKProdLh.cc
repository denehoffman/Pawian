#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKProdLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh(boost::shared_ptr<const JpsiGamKsKlKKEventList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsJpsiGamKsKlKKLh(theEvtList)
  ,_eta2225Hyp(true)
  ,_f02020Hyp(true)
  ,_nFitParams(0)    
{

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

  _ampVec.push_back(paramEnumGamKsKlKK::etacGamma);
  _massVec.push_back(paramEnumGamKsKlKK::etac);
  _nFitParams+=4;

  if(_eta2225Hyp){
    _ampVec.push_back(paramEnumGamKsKlKK::eta2225Gamma); 
    _massVec.push_back(paramEnumGamKsKlKK::eta2225);
    _nFitParams+=4;
  }

  if(_f02020Hyp){
    _ampVec.push_back(paramEnumGamKsKlKK::f02020Gamma);
    _ampVec.push_back(paramEnumGamKsKlKK::f02020ToPhiPhi); 
    _massVec.push_back(paramEnumGamKsKlKK::f02020);
    _nFitParams+=10;
  }
//   _ampVec.push_back(paramEnumGamKsKlKK::f22010Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22300Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22340Gamma);
 
  
//   std::vector<unsigned int>::iterator ampIt;
//   for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
//     std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitparamsGamKsKlKK.jpclsVec(*ampIt);
// //     _nFitParams+=2*JPCLSs.size();
//   }

 //  _nFitParams+=4; //2 masses+ 2 widths
//   _nFitParams= int(_ampVec.size())*2+int(_massVec.size())*2;
  
}

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh( boost::shared_ptr<AbsJpsiGamKsKlKKLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsJpsiGamKsKlKKLh(theLhPtr->getEventList())
  ,_eta2225Hyp(true)
  ,_nFitParams(0) 
{

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

  _ampVec.push_back(paramEnumGamKsKlKK::etacGamma);
  _massVec.push_back(paramEnumGamKsKlKK::etac);
  _nFitParams=4;

  if(_eta2225Hyp){
    _ampVec.push_back(paramEnumGamKsKlKK::eta2225Gamma);
    _massVec.push_back(paramEnumGamKsKlKK::eta2225);
    _nFitParams+=4;
  }

  if(_f02020Hyp){
    _ampVec.push_back(paramEnumGamKsKlKK::f02020Gamma);
    _ampVec.push_back(paramEnumGamKsKlKK::f02020ToPhiPhi); 
    _massVec.push_back(paramEnumGamKsKlKK::f02020);
    _nFitParams+=10;
  }
//   _ampVec.push_back(paramEnumGamKsKlKK::f22010Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22300Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22340Gamma);
  

  

//   std::vector<unsigned int>::iterator ampIt;
//   for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
//     std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitparamsGamKsKlKK.jpclsVec(*ampIt);
//     _nFitParams+=2*JPCLSs.size();
//   }

}

JpsiGamKsKlKKProdLh::~JpsiGamKsKlKKProdLh()
{;
}



double JpsiGamKsKlKKProdLh::calcEvtIntensity(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, const paramGamKsKlKK& theParamVal){

  double result=0.;

  Vector4<double> fvKsKlKpKm( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );

//   complex<double> JmpGmp=etacGammaCoherentAmp(1, 0, 1, theParamVal, theData);
//   complex<double> JmpGmm=etacGammaCoherentAmp(1, 0, -1, theParamVal, theData);
//   complex<double> JmmGmp=etacGammaCoherentAmp(-1, 0, 1, theParamVal, theData);
//   complex<double> JmmGmm=etacGammaCoherentAmp(-1, 0, -1, theParamVal, theData);

//   if (_eta2225Hyp){
//     JmpGmp+=eta2225GammaCoherentAmp(1, 0, 1, theParamVal, theData);
//     JmpGmm+=eta2225GammaCoherentAmp(1, 0, -1, theParamVal, theData);
//     JmmGmp+=eta2225GammaCoherentAmp(-1, 0, 1, theParamVal, theData);
//     JmmGmm+=eta2225GammaCoherentAmp(-1, 0, -1, theParamVal, theData);
//   }

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEtacGam=theParamVal.PsiToEtacGamma;
  pair<double, double> BwEtacParam=theParamVal.BwEtac;
  complex<double> JmpGmp=etaGammaAmp(1, 0, 1, theData, PsiToEtacGam, BwEtacParam);
  complex<double> JmpGmm=etaGammaAmp(1, 0, -1, theData, PsiToEtacGam, BwEtacParam);
  complex<double> JmmGmp=etaGammaAmp(-1, 0, 1, theData, PsiToEtacGam, BwEtacParam);
  complex<double> JmmGmm=etaGammaAmp(-1, 0, -1, theData, PsiToEtacGam, BwEtacParam); 

  if (_eta2225Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEta2225Gam=theParamVal.PsiToEta2225Gamma;
    pair<double, double> BwEta2225Param=theParamVal.BwEta2225;
    JmpGmp+=etaGammaAmp(1, 0, 1, theData, PsiToEta2225Gam, BwEta2225Param );
    JmpGmm+=etaGammaAmp(1, 0, -1, theData, PsiToEta2225Gam, BwEta2225Param);
    JmmGmp+=etaGammaAmp(-1, 0, 1, theData, PsiToEta2225Gam, BwEta2225Param);
    JmmGmm+=etaGammaAmp(-1, 0, -1, theData, PsiToEta2225Gam, BwEta2225Param);
  }

  if (_f02020Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiTof02020Gam=theParamVal.PsiToF02020Gamma;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > F02020ToPhiPhi=theParamVal.F02020ToPhiPhi;
    pair<double, double> BwF02020Param=theParamVal.BwF02020;
    JmpGmp+=f0GammaAmp(1, 1, theData, PsiTof02020Gam, F02020ToPhiPhi, BwF02020Param );
    JmpGmm+=f0GammaAmp(1, -1, theData, PsiTof02020Gam, F02020ToPhiPhi, BwF02020Param );
    JmmGmp+=f0GammaAmp(-1, 1, theData, PsiTof02020Gam, F02020ToPhiPhi, BwF02020Param );
    JmmGmm+=f0GammaAmp(-1, -1, theData, PsiTof02020Gam, F02020ToPhiPhi, BwF02020Param );
  } 

  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm); 
  return result;  
}

complex<double> JpsiGamKsKlKKProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}


complex<double> JpsiGamKsKlKKProdLh::etaGammaAmp(Spin Minit, Spin Metac, Spin Mgamma, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampParam, pair<double, double>& massParam){

   complex<double> result(0.,0.);

   Vector4<double> fv2Phi( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=ampParam.begin(); itPsi!=ampParam.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second.first;
     double thePsiPhi=itPsi->second.second;
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = Metac-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(0, Metac, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->df_Psi[PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }

   result*=BreitWigner( fv2Phi, massParam.first, massParam.second)*etaToPhiPhiTo4KAmp(theData);

return result;
}

complex<double> JpsiGamKsKlKKProdLh::f0GammaAmp(Spin Minit, Spin Mgamma, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampf0Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampf0Dec, pair<double, double>& massParam){

   complex<double> result(0.,0.);
   Spin f0Spin=0;
   Vector4<double> fv2Phi( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=ampf0Prod.begin(); itPsi!=ampf0Prod.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second.first;
     double thePsiPhi=itPsi->second.second;
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = f0Spin-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(0, f0Spin, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->df_Psi[PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }

   result*=BreitWigner( fv2Phi, massParam.first, massParam.second)*f0ToPhiPhiTo4KAmp(theData, ampf0Dec);

return result;
}


complex<double> JpsiGamKsKlKKProdLh::f22010GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  
   complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEtacGam=theParamVal.PsiToF22010Gamma;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToEtacGam.begin(); itPsi!=PsiToEtacGam.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second.first;
     double thePsiPhi=itPsi->second.second;
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = Metac-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(2, Metac, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->df_Psi[PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }
   Vector4<double> fv( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   result*=BreitWigner( fv, theParamVal.BwF22010.first, theParamVal.BwF22010.second);
   //result*=BreitWigner( fv, 2.011   ,0.202 );
   return result;
}   



complex<double> JpsiGamKsKlKKProdLh::f22300GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  
   complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEtacGam=theParamVal.PsiToF22300Gamma;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToEtacGam.begin(); itPsi!=PsiToEtacGam.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second.first;
     double thePsiPhi=itPsi->second.second;
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = Metac-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(2, Metac, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->df_Psi[PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }
   Vector4<double> fv( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   result*=BreitWigner( fv, theParamVal.BwF22300.first, theParamVal.BwF22300.second);
   //result*=BreitWigner( fv, 2.297,   0.149 );
   return result;
}   






complex<double> JpsiGamKsKlKKProdLh::f22340GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  
   complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEtacGam=theParamVal.PsiToF22340Gamma;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToEtacGam.begin(); itPsi!=PsiToEtacGam.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second.first;
     double thePsiPhi=itPsi->second.second;
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
     Spin lambda = Metac-Mgamma;
     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(2, Metac, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->df_Psi[PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }
   Vector4<double> fv( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   result*=BreitWigner( fv, theParamVal.BwF22340.first, theParamVal.BwF22340.second);
   //result*=BreitWigner( fv, 2.339   ,0.319 );
   return result;
}   


complex<double> JpsiGamKsKlKKProdLh::etaToPhiPhiTo4KAmp(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){

  complex<double> result(0.,0.);

    for (Spin lambdaPhi=-1; lambdaPhi<1; lambdaPhi++){
       result+=Clebsch(1, lambdaPhi, 1, -lambdaPhi, 1, 0)*3.*conj(theData->Df_KsKl[1][lambdaPhi][0])*conj(theData->Df_KpKm[1][lambdaPhi][0]); //Clebsch(1, lambdaPhi, 1, -lambdaPhi, 1, 0): jj-coupling mother->phi phi; 3=sqrt(2L+1)*sqrt(2L+1) Cls=1


      }

    return result;
}

complex<double> JpsiGamKsKlKKProdLh::f0ToPhiPhiTo4KAmp(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampf0Dec){
  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itf0;
   for ( itf0=ampf0Dec.begin(); itf0!=ampf0Dec.end(); ++itf0){
     boost::shared_ptr<const JPCLS> f0ToPhiPhiState=itf0->first;
     double theMag=itf0->second.first;
     double thePhi=itf0->second.second;
     complex<double> expiphi(cos(thePhi), sin(thePhi));

     complex<double> tmp2PhiDecAmp(0.,0.);
     for (Spin lambdaPhi1=-1; lambdaPhi1<1; lambdaPhi1++){
       Spin lambdaPhi2=lambdaPhi1;
       tmp2PhiDecAmp+=Clebsch(1, lambdaPhi1, 1, -lambdaPhi2, f0ToPhiPhiState->S, lambdaPhi1-lambdaPhi2)
	 *3.*conj(theData->Df_KsKl[1][lambdaPhi1][0])*conj(theData->Df_KpKm[1][lambdaPhi2][0]); //3=sqrt(2L+1)*sqrt(2L+1) Cls=1    
     }
   
     result+= theMag*expiphi*sqrt(2*f0ToPhiPhiState->L+1)*tmp2PhiDecAmp;
   }
   return result;
}


void JpsiGamKsKlKKProdLh::setMnUsrParams(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal){

//   _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::etacGamma);
//   _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::eta2225Gamma);
//   _fitparamsGamKsKlKK.setMnUsrParamsMass(upar,startVal,errVal,paramEnumGamKsKlKK::etac);
//   _fitparamsGamKsKlKK.setMnUsrParamsMass(upar,startVal,errVal,paramEnumGamKsKlKK::eta2225);

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitparamsGamKsKlKK.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    _fitparamsGamKsKlKK.setMnUsrParamsMass(upar, startVal, errVal, (*itMasses) );  
  }

//   _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::f22010Gamma);
//   _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::f22300Gamma);
//   _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::f22340Gamma);
}



int JpsiGamKsKlKKProdLh::setFitParamVal(paramGamKsKlKK& theParamVal, const std::vector<double>& par){  

  
  if (par.size()!= nFitParams() ) {
    Alert  << "size of parameters wrong!!! par.size()=" << par.size() << 
    "; it should be " << nFitParams() ;  // << endmsg;
    exit(1);
  }
  
  int counter=0;

//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   //Psi(2S) ->Chi_c1 gamma amplitude params
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::etacGamma);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::eta2225Gamma);
//    counter=_fitparamsGamKsKlKK.setFitParamValMass(theParamVal, par, counter,  paramEnumGamKsKlKK::etac);
//   counter=_fitparamsGamKsKlKK.setFitParamValMass(theParamVal, par, counter, paramEnumGamKsKlKK::eta2225);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::f22010Gamma);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::f22300Gamma);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::f22340Gamma);


  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    counter=_fitparamsGamKsKlKK.setFitParamValMass(theParamVal, par, counter, (*itMasses) ); 
  }


  return counter;
}

unsigned int  JpsiGamKsKlKKProdLh::nFitParams(){
  return _nFitParams;
}

void JpsiGamKsKlKKProdLh::print(std::ostream& os) const{
  os << "JpsiGamKsKlKKProdLh::print\n";
}

void JpsiGamKsKlKKProdLh::printCurrentFitResult(paramGamKsKlKK& theParamVal){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitparamsGamKsKlKK.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      Info << "(debug)" <<(*itJPCLS)->name()<< paramEnumGamKsKlKK::name(*itAmps);// ;  // << endmsg;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitparamsGamKsKlKK.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      Info <<"(debug)\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second;//  ;  // << endmsg;
    }  
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    DebugMsg<< paramEnumGamKsKlKK::name(*itMasses) << endmsg;
    std::pair<double, double> tmpParam=_fitparamsGamKsKlKK.massPair(theParamVal, *itMasses);
    DebugMsg <<"\t Mass:" << tmpParam.first <<"\t Width:" << tmpParam.second  << endmsg;
  }

}

void JpsiGamKsKlKKProdLh::dumpCurrentResult(std::ostream& os, paramGamKsKlKK& theParamVal, std::string& suffix){

  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" ;  // << endmsg;
    return;
  }

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitparamsGamKsKlKK.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitparamsGamKsKlKK.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];

      std::string tmpStringDec=(*itJPCLS)->name()+paramEnumGamKsKlKK::name(*itAmps)+suffix;
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }  
  }

}
