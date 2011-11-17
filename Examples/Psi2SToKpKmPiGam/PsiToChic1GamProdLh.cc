#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/PsiToChic1GamProdLh.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToChic1GamProdLh::PsiToChic1GamProdLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  AbsLh(theEvtList)
{
}

PsiToChic1GamProdLh::PsiToChic1GamProdLh( boost::shared_ptr<AbsLh> theLhPtr ) :
  AbsLh(theLhPtr->getEventList())
{
}

PsiToChic1GamProdLh::~PsiToChic1GamProdLh()
{;
}



double PsiToChic1GamProdLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){

  double result=0.;
  
   for (Spin PsiM=-1; PsiM<=1; PsiM=PsiM+2){
    for (Spin ChiM=-1; ChiM<=1; ChiM++){
      for (Spin GamM=-1; GamM<=1; GamM=GamM+2){

	result+=norm(calcCoherentAmp(PsiM, ChiM, GamM, theParamVal, theData));
      }
    }
  }
  return result;  
}

complex<double> PsiToChic1GamProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}

complex<double> PsiToChic1GamProdLh::calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, fitParams& theParamVal, EvtData* theData){
  
  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& thePsiMags=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::ChiGam];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& thePsiPhis=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::ChiGam]; 

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsiMag;
  for ( itPsiMag=thePsiMags.begin(); itPsiMag!=thePsiMags.end(); ++itPsiMag){
    
    boost::shared_ptr<const JPCLS> PsiState=itPsiMag->first;
    double thePsiMag=itPsiMag->second;
    double thePsiPhi=thePsiPhis[PsiState];
    complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
    
    Spin lamChiGam=MChi-MGamma;
    if (fabs(lamChiGam)>PsiState->J || fabs(lamChiGam)>PsiState->S) continue;
    
    complex<double> PsiDecAmp=thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
      *Clebsch(PsiState->L,0,PsiState->S,lamChiGam,PsiState->J,lamChiGam)
      *Clebsch(1, MChi, 1, -MGamma, PsiState->S,lamChiGam)
      *conj(theData->WignerDs[enumChic1ToKpKmPiGamData::dfPsi][PsiState->J][Minit][lamChiGam]);
    
    result+=PsiDecAmp;
  }      

  return result; 
}

void PsiToChic1GamProdLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  FitParamsChic1ToKpKmPiGam fitParamChic1;


  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;
  theAmpMap[paramEnumChic1ToKpKmPiGam::ChiGam] = fitParamChic1.jpclsVec(paramEnumChic1ToKpKmPiGam::ChiGam);
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > >::iterator itAmpMap;

  for (itAmpMap=theAmpMap.begin(); itAmpMap!=theAmpMap.end(); ++itAmpMap){

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > valMagMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > errMagMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > valPhiMap;
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > errPhiMap;
    
    std::vector< boost::shared_ptr<const JPCLS> >::iterator itAmp;
    for (itAmp=itAmpMap->second.begin(); itAmp!=itAmpMap->second.end(); ++itAmp){
        valMagMap[(*itAmp)]=0.5;
        errMagMap[(*itAmp)]=0.2; 
        valPhiMap[(*itAmp)]=0.;
        errPhiMap[(*itAmp)]=0.2;      
      }

    fitVal.Mags[itAmpMap->first]=valMagMap;
    fitVal.Phis[itAmpMap->first]=valPhiMap;  
    fitErr.Mags[itAmpMap->first]=errMagMap;
    fitErr.Phis[itAmpMap->first]=errPhiMap;  
  }

}

void PsiToChic1GamProdLh::print(std::ostream& os) const{
  os << "PsiToChic1GamProdLh::print\n";
}

