#include <fstream>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKStreamFitParams.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKStreamFitParams::JpsiGamKsKlKKStreamFitParams(std::string& filePath) :
  AbsFitParamStreamer(filePath),
  _fitparamsGamKsKlKK()
{
  fillParamMap();

}

JpsiGamKsKlKKStreamFitParams::~JpsiGamKsKlKKStreamFitParams(){;}

void JpsiGamKsKlKKStreamFitParams::fillParamMap(){
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;


  for (unsigned int ui=paramEnumGamKsKlKK::etacGamma; ui<paramEnumGamKsKlKK::nAmps; ui++){
    
    std::string theAmpString=paramEnumGamKsKlKK::name(ui);
    std::string theSuffix=theAmpString+"Val";
    std::vector< boost::shared_ptr<const JPCLS> >  theJPCLSs=_fitparamsGamKsKlKK.jpclsVec(ui);

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& theAmpMapVal=_fitparamsGamKsKlKK.ampMap(_paramVal, ui);


    fillParamMapAmps(theJPCLSs, theSuffix, theAmpMapVal);

    theSuffix=theAmpString+"Err";
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& theAmpMapErr=_fitparamsGamKsKlKK.ampMap(_paramErr, ui);
    fillParamMapAmps(theJPCLSs, theSuffix, theAmpMapErr);
  }


  for (unsigned int ui=paramEnumGamKsKlKK::etac; ui<paramEnumGamKsKlKK::nMasses; ui++){
    std::string theMassString=paramEnumGamKsKlKK::name(ui);
    std::string theSuffix=theMassString+"Val";

    pair<double, double>& theValPairs=_fitparamsGamKsKlKK.massPair(_paramVal, ui);
    fillParamMapMass(theSuffix, theValPairs);

    theSuffix=theMassString+"Err";
    pair<double, double>& theErrPairs=_fitparamsGamKsKlKK.massPair(_paramErr, ui);
    fillParamMapMass(theSuffix, theErrPairs);

  }

}
