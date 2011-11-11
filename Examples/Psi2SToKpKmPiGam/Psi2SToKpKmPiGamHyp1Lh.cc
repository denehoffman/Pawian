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



void Psi2SToKpKmPiGamHyp1Lh::print(std::ostream& os) const{
  os << "Psi2SToKpKmPiGamHyp1Lh::print\n";
}

