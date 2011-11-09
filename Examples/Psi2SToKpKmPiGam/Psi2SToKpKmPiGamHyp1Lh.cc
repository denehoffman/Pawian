#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamHyp1Lh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamHyp1Lh::Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  Psi2SToKpKmPiGamBaseLh(theEvtList,  hypMap)
{
}

Psi2SToKpKmPiGamHyp1Lh::Psi2SToKpKmPiGamHyp1Lh( boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Psi2SToKpKmPiGamBaseLh(theLhPtr->getEventList(), hypMap)
{
}

Psi2SToKpKmPiGamHyp1Lh::~Psi2SToKpKmPiGamHyp1Lh()
{;
}


complex<double> Psi2SToKpKmPiGamHyp1Lh::calcDecAmp(complex<double>& inAmp,Spin lamChi, fitParams& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){

  
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

