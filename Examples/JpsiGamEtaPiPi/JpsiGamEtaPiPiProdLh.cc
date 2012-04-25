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
  ,_etaHyp(false)
  ,_f1Hyp(false)
  ,_usePhasespace(false)
{
  initializeHypothesisMap( hypMap);
 
}

JpsiGamEtaPiPiProdLh::JpsiGamEtaPiPiProdLh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsLh(theLhPtr->getEventList())
  ,_etaHyp(false)
  ,_f1Hyp(false)
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
  
  if(_etaHyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiToEtacGamPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma];
    
    calcEtaGammaAmp( theData, PsiToEtacGamMag, PsiToEtacGamPhi, JmpGmp, JmpGmm, JmmGmp, JmmGmm);
    
  }

  if(_f1Hyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof1GamMag=theParamVal.Mags[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > PsiTof1GamPhi=theParamVal.Phis[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma];
//     JmpGmp+=etaGammaAmp(1, 0, 1, theData, PsiToEtaGamMag, PsiToEtaGamPhi );
//     JmpGmm+=etaGammaAmp(1, 0, -1, theData,  PsiToEtaGamMag, PsiToEtaGamPhi  );
//     JmmGmp+=etaGammaAmp(-1, 0, 1, theData,  PsiToEtaGamMag, PsiToEtaGamPhi );
//     JmmGmm+=etaGammaAmp(-1, 0, -1, theData,  PsiToEtaGamMag, PsiToEtaGamPhi  );
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

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;
   for ( itPsi=PsiToXGamMag.begin(); itPsi!=PsiToXGamMag.end(); ++itPsi){
     boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
     double thePsiMag=itPsi->second;
     double thePsiPhi=PsiToXGamPhi[PsiState];
     complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));

     Spin lambda = lamX-lamGamma;

     if( fabs(lambda)>PsiState->J || fabs(lambda)>PsiState->S) continue;
     
     complex<double> amp = thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
       *Clebsch(PsiState->L, 0, PsiState->S, lambda, PsiState->J, lambda)
       *Clebsch(jX, lamX, 1, -lamGamma, PsiState->S, lambda  )
       *conj( theData->WignerDs[enumJpsiGamEtaPiPiData::Df_Psi][PsiState->J][Minit][lambda]  );
     
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





void JpsiGamEtaPiPiProdLh::print(std::ostream& os) const{
  os << "JpsiGamEtaPiPiProdLh::print\n";
}


void JpsiGamEtaPiPiProdLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){
  JpsiGamEtaPiPiFitParams theFitParams;
  
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;
  
  if(_etaHyp){
    theAmpMap[paramEnumJpsiGamEtaPiPi::PsiToEtaGamma] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::PsiToEtaGamma);
  }
  
  if(_f1Hyp){
    theAmpMap[paramEnumJpsiGamEtaPiPi::PsiToF1Gamma] = theFitParams.jpclsVec(paramEnumJpsiGamEtaPiPi::PsiToF1Gamma);
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

  //fill other params  
  if(_usePhasespace){
    fitVal.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace]=0.2;
    fitErr.otherParams[paramEnumJpsiGamEtaPiPi::phaseSpace]=0.4;
  }

}
  



bool 
JpsiGamEtaPiPiProdLh::initializeHypothesisMap( const std::map<const std::string, bool>& hypMap   ){
  
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("etaHyp");
  
  if (iter !=hypMap.end()){
    _etaHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _etaHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis etaHyp not set!!!" <<endmsg;
  
  
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








