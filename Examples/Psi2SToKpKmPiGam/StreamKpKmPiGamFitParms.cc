#include <fstream>

#include "Examples/Psi2SToKpKmPiGam/StreamKpKmPiGamFitParms.hh"
#include "ErrLogger/ErrLogger.hh"

StreamKpKmPiGamFitParms::StreamKpKmPiGamFitParms(std::string& filePath) :
  AbsFitParamStreamer(filePath),
  _fitParamsKpKmPiGam()
{
  fillParamMap();

}

StreamKpKmPiGamFitParms::~StreamKpKmPiGamFitParms(){;}

void StreamKpKmPiGamFitParms::fillParamMap(){
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;


  for (unsigned int ui=paramEnumKpKmPiGam::ChiGam; ui<paramEnumKpKmPiGam::nAmps; ui++){

    std::string theAmpString=paramEnumKpKmPiGam::name(ui);
    std::string theSuffix=theAmpString+"Val";
    std::vector< boost::shared_ptr<const JPCLS> >  theJPCLSs=_fitParamsKpKmPiGam.jpclsVec(ui);

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& theAmpMapVal=_fitParamsKpKmPiGam.ampMap(_paramVal, ui);


    fillParamMapAmps(theJPCLSs, theSuffix, theAmpMapVal);

    theSuffix=theAmpString+"Err";
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& theAmpMapErr=_fitParamsKpKmPiGam.ampMap(_paramErr, ui);
    fillParamMapAmps(theJPCLSs, theSuffix, theAmpMapErr);
  }


  for (unsigned int ui=paramEnumKpKmPiGam::K890; ui<paramEnumKpKmPiGam::nMasses; ui++){
    std::string theMassString=paramEnumKpKmPiGam::name(ui);
    std::string theSuffix=theMassString+"Val";

    pair<double, double>& theValPairs=_fitParamsKpKmPiGam.massPair(_paramVal, ui);
    fillParamMapMass(theSuffix, theValPairs);

    theSuffix=theMassString+"Err";
    pair<double, double>& theErrPairs=_fitParamsKpKmPiGam.massPair(_paramErr, ui);
    fillParamMapMass(theSuffix, theErrPairs);

  }


  std::string theKey="FlatMa980Val";
  fillParamFlatte(theKey, _paramVal.FlatMa980);
  theKey="FlatMa980Err";
  fillParamFlatte(theKey, _paramErr.FlatMa980);

  theKey="FlatgKKVal";
  fillParamFlatte(theKey, _paramVal.FlatgKK);
  theKey="FlatgKKErr";
  fillParamFlatte(theKey, _paramErr.FlatgKK);

  theKey="FlatgEtaPiVal";
  fillParamFlatte(theKey, _paramVal.FlatgEtaPi);
  theKey="FlatgEtaPiErr";
  fillParamFlatte(theKey, _paramErr.FlatgEtaPi);

  theKey="phaseSpaceVal";
  fillParamFlatte(theKey, _paramVal.phaseSpace);
  theKey="phaseSpaceErr";
  fillParamFlatte(theKey, _paramErr.phaseSpace);

}
