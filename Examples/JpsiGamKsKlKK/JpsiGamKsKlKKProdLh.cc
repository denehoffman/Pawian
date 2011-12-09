#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKProdLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh(boost::shared_ptr<const JpsiGamKsKlKKEventList> theEvtList) :
  AbsJpsiGamKsKlKKLh(theEvtList)
  ,_nFitParams(0)    
{
  _ampVec.push_back(paramEnumGamKsKlKK::etacGamma);
  _ampVec.push_back(paramEnumGamKsKlKK::eta2225Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22010Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22300Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22340Gamma);
 
  
  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitparamsGamKsKlKK.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  _nFitParams+=4; //2 masses+ 2 widths
  
}

JpsiGamKsKlKKProdLh::JpsiGamKsKlKKProdLh( boost::shared_ptr<AbsJpsiGamKsKlKKLh> theLhPtr ) :
  AbsJpsiGamKsKlKKLh(theLhPtr->getEventList())
  ,_nFitParams(0) 
{
  _ampVec.push_back(paramEnumGamKsKlKK::etacGamma);
  _ampVec.push_back(paramEnumGamKsKlKK::eta2225Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22010Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22300Gamma);
//   _ampVec.push_back(paramEnumGamKsKlKK::f22340Gamma);
  

  

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitparamsGamKsKlKK.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }
  _nFitParams+=4; //2 masses+ 2 widths
}

JpsiGamKsKlKKProdLh::~JpsiGamKsKlKKProdLh()
{;
}



double JpsiGamKsKlKKProdLh::calcEvtIntensity(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, const paramGamKsKlKK& theParamVal){

  double result=0.;

  complex<double> JmpGmp=etacGammaCoherentAmp(1, 0, 1, theParamVal, theData)+eta2225GammaCoherentAmp(1, 0, 1, theParamVal, theData);
  complex<double> JmpGmm=etacGammaCoherentAmp(1, 0, -1, theParamVal, theData)+eta2225GammaCoherentAmp(1, 0, -1, theParamVal, theData);
  complex<double> JmmGmp=etacGammaCoherentAmp(-1, 0, 1, theParamVal, theData)+eta2225GammaCoherentAmp(-1, 0, 1, theParamVal, theData);
  complex<double> JmmGmm=etacGammaCoherentAmp(-1, 0, -1, theParamVal, theData)+eta2225GammaCoherentAmp(-1, 0, -1, theParamVal, theData);

  
//   for (Spin PsiM=-1; PsiM<=1; PsiM=PsiM+2){
//     for (Spin GamM=-1; GamM<=1; GamM=GamM+2){
//       complex<double> tmp(0.,0.);
//       //tmp+=(eta2225GammaCoherentAmp(PsiM, 0, GamM, theParamVal, theData));
//       tmp+=(etacGammaCoherentAmp(PsiM, 0, GamM, theParamVal, theData));
      
//       for(Spin f2M=-2; f2M<=2; f2M++){
// 	//      tmp+=(f22340GammaCoherentAmp(PsiM, f2M, GamM, theParamVal, theData));
// 	//	tmp+=(f22300GammaCoherentAmp(PsiM, f2M, GamM, theParamVal, theData));
// 	tmp+=(f22010GammaCoherentAmp(PsiM, f2M, GamM, theParamVal, theData));
//       }
      
//       result+=tmp;
//     }
//   }
  result=norm(JmpGmp)+norm(JmpGmm)+norm(JmmGmp)+norm(JmmGmm); 
  return result;  
}

complex<double> JpsiGamKsKlKKProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}





complex<double> JpsiGamKsKlKKProdLh::etacGammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  
   complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEtacGam=theParamVal.PsiToEtacGamma;
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
       *Clebsch(0, Metac, 1, -Mgamma, PsiState->S, lambda  )
       *conj( theData->df_Psi[PsiState->J][Minit][lambda]  );
     
     result+= amp;
   }
   Vector4<double> fv( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   
   //cout << "Parm etac " << theParamVal.BwEtac.first << " " << theParamVal.BwEtac.second << endl;

   result*=BreitWigner( fv, theParamVal.BwEtac.first, theParamVal.BwEtac.second);
   //result*=BreitWigner( fv,2.98 , 0.027);
return result;
}   

complex<double> JpsiGamKsKlKKProdLh::eta2225GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){
  
   complex<double> result(0.,0.);
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToEta2225Gam=theParamVal.PsiToEta2225Gamma;
   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToEta2225Gam.begin(); itPsi!=PsiToEta2225Gam.end(); ++itPsi){
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
   Vector4<double> fv( theData->V4_KsKlKpKm_HeliPsi.E(), theData->V4_KsKlKpKm_HeliPsi.X(), theData->V4_KsKlKpKm_HeliPsi.Y(), theData->V4_KsKlKpKm_HeliPsi.Z() );
   result*=BreitWigner( fv, theParamVal.BwEta2225.first, theParamVal.BwEta2225.second);
   //result*=BreitWigner( fv, 2.226, 0.185);
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



void JpsiGamKsKlKKProdLh::setMnUsrParams(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal){

  _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::etacGamma);
  _fitparamsGamKsKlKK.setMnUsrParamsDec(upar,startVal,errVal,paramEnumGamKsKlKK::eta2225Gamma);
  _fitparamsGamKsKlKK.setMnUsrParamsMass(upar,startVal,errVal,paramEnumGamKsKlKK::etac);
  _fitparamsGamKsKlKK.setMnUsrParamsMass(upar,startVal,errVal,paramEnumGamKsKlKK::eta2225);
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

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::etacGamma);
  counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::eta2225Gamma);
   counter=_fitparamsGamKsKlKK.setFitParamValMass(theParamVal, par, counter,  paramEnumGamKsKlKK::etac);
  counter=_fitparamsGamKsKlKK.setFitParamValMass(theParamVal, par, counter, paramEnumGamKsKlKK::eta2225);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::f22010Gamma);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::f22300Gamma);
//   counter=_fitparamsGamKsKlKK.setFitParamValDec(theParamVal, par, counter, paramEnumGamKsKlKK::f22340Gamma);



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
