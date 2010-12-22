#include <fstream>

#include "Examples/Psi2STo2K2PiGam/Stream2K2PiGamFitParms.hh"
#include "ErrLogger/ErrLogger.hh"

Stream2K2PiGamFitParms::Stream2K2PiGamFitParms(std::string& filePath) :
  AbsFitParamStreamer(filePath),
  _fitParams2K2PiGam()
{
  fillParamMap();

}

Stream2K2PiGamFitParms::~Stream2K2PiGamFitParms(){;}

void Stream2K2PiGamFitParms::fillParamMap(){
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;


  for (unsigned int ui=paramEnum2K2PiGam::ChiGam; ui<paramEnum2K2PiGam::nAmps; ui++){

    std::string theAmpString=paramEnum2K2PiGam::name(ui);
    std::string theSuffix=theAmpString+"Val";
    std::vector< boost::shared_ptr<const JPCLS> >  theJPCLSs=_fitParams2K2PiGam.jpclsVec(ui);

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& theAmpMapVal=_fitParams2K2PiGam.ampMap(_paramVal, ui);


    fillParamMapAmps(theJPCLSs, theSuffix, theAmpMapVal);

    theSuffix=theAmpString+"Err";
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& theAmpMapErr=_fitParams2K2PiGam.ampMap(_paramErr, ui);
    fillParamMapAmps(theJPCLSs, theSuffix, theAmpMapErr);
  }


  for (unsigned int ui=paramEnum2K2PiGam::K892; ui<paramEnum2K2PiGam::nMasses; ui++){
    std::string theMassString=paramEnum2K2PiGam::name(ui);
    std::string theSuffix=theMassString+"Val";

    pair<double, double>& theValPairs=_fitParams2K2PiGam.massPair(_paramVal, ui);
    fillParamMapMass(theSuffix, theValPairs);

    theSuffix=theMassString+"Err";
    pair<double, double>& theErrPairs=_fitParams2K2PiGam.massPair(_paramErr, ui);
    fillParamMapMass(theSuffix, theErrPairs);

  }


  std::string theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980)+"Val";
  fillParamFlatte(theKey, _paramVal.Flatf980);
  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980)+"Err";
  fillParamFlatte(theKey, _paramErr.Flatf980);

  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gKK)+"Val";
  fillParamFlatte(theKey, _paramVal.Flatf980gKK);
  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gKK)+"Err";
  fillParamFlatte(theKey, _paramErr.Flatf980gKK);

  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gpipi)+"Val";
  fillParamFlatte(theKey, _paramVal.Flatf980gPiPi);
  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gpipi)+"Err";
  fillParamFlatte(theKey, _paramErr.Flatf980gPiPi);

  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace)+"Val";
  fillParamFlatte(theKey, _paramVal.phaseSpace);
  theKey=paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace)+"Err";
  fillParamFlatte(theKey, _paramErr.phaseSpace);

}
