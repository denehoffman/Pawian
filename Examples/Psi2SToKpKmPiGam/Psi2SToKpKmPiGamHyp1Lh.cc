#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamHyp1Lh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEvtList.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamHyp1Lh::Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<const EvtDataBaseList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  Psi2SToKpKmPiGamBaseLh(theEvtList,  hypMap)
{
}

Psi2SToKpKmPiGamHyp1Lh::Psi2SToKpKmPiGamHyp1Lh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Psi2SToKpKmPiGamBaseLh(theLhPtr->getEventList(), hypMap)
{
}

Psi2SToKpKmPiGamHyp1Lh::~Psi2SToKpKmPiGamHyp1Lh()
{;
}


complex<double> Psi2SToKpKmPiGamHyp1Lh::calcDecAmp(complex<double>& inAmp,Spin lamChi, fitParams& theParamVal, EvtData* theData){

  
   complex<double> result=Psi2SToKpKmPiGamBaseLh::calcDecAmp(inAmp, lamChi, theParamVal, theData);

   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToKappaKMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::KappaK];
   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToKappaKPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::KappaK];

  double KappaMass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::Kappa];
  double KappaWidth=theParamVal.Widths[paramEnumChic1ToKpKmPiGam::Kappa];

  //Chi to a0(980) pi0 and a0(980)->K+K-
  result+=K0_Amp(theData, ChiToKappaKMag,  ChiToKappaKPhi, KappaMass, KappaWidth, lamChi);
// //   result+=KappaPole_Amp(theData, ChiToKappaK, KappaMass, KappaWidth, lamChi);
  return result;
}



void Psi2SToKpKmPiGamHyp1Lh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){

  Psi2SToKpKmPiGamBaseLh::getDefaultParams(fitVal, fitErr);

  FitParamsChic1ToKpKmPiGam fitParamChic1;
  std::map<int, std::vector< boost::shared_ptr<const JPCLS> > > theAmpMap;

  theAmpMap[paramEnumChic1ToKpKmPiGam::KappaK]=fitParamChic1.jpclsVec(paramEnumChic1ToKpKmPiGam::KappaK);

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

  //fill masses and widths
  fitVal.Masses[paramEnumChic1ToKpKmPiGam::Kappa]=0.6;    
  fitErr.Masses[paramEnumChic1ToKpKmPiGam::Kappa]=0.1;
  fitVal.Widths[paramEnumChic1ToKpKmPiGam::Kappa]=0.6;
  fitErr.Widths[paramEnumChic1ToKpKmPiGam::Kappa]=0.1;

  //fill gFactors

  //fill others
}



void Psi2SToKpKmPiGamHyp1Lh::print(std::ostream& os) const{
  os << "Psi2SToKpKmPiGamHyp1Lh::print\n";
}

