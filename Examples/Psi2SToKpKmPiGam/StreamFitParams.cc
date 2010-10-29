#include <fstream>

#include "Examples/Psi2SToKpKmPiGam/StreamFitParams.hh"
#include "ErrLogger/ErrLogger.hh"

StreamFitParams::StreamFitParams(std::string& filePath, boost::shared_ptr<const Psi2SToKpKmPiGamStates> theStates)
{

//   std::ifstream ifs(filePath.c_str(), boost::shared_ptr<const Psi2SToKpKmPiGamStates>);
  std::ifstream ifs(filePath.c_str());
  if(!ifs.good()) 
    {
      Alert << "Error: could not parse " << filePath << endmsg;
      exit(1);
    }


  string strTmp, strState;
  bool bError=false;
  double dMag=0.0, dPhi=0.0;
  
  while(!ifs.eof())
    {
      ifs >> strState;
      if (ifs.eof()) break;
      ifs >> strTmp;
      dMag = atof(strTmp.c_str());
      if (ifs.eof()) break;
      ifs >> strTmp;
      dPhi = atof(strTmp.c_str());
      if (ifs.eof()) break;
      if (!bError) 
	{
	  DebugMsg << "Start Parameter: " << strState << " Mag: " << dMag << " Phi: " << dPhi << " set.\n" << endmsg;
	  StringParameterMap(strState,dMag,dPhi);
	}
    }

  fillParamMap(theStates);

}

void StreamFitParams::StringParameterMap(const string &theName, const double firstVal, const double secondVal)
{
  std::pair <double,double> thePair=make_pair(firstVal, secondVal);
  _stringPairMap[theName]=thePair;
}

void StreamFitParams::fillParamMap(boost::shared_ptr<const Psi2SToKpKmPiGamStates> psi2SToKpKmPiGamStatesPtr){
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=psi2SToKpKmPiGamStatesPtr->PsiToChiGamStates();
 
  std::string theSuffix("ChiGamVal");
  fillParamMapAmps(JPCLSPsiToChiGam, theSuffix, _paramVal.PsiToChiGam);
  theSuffix="ChiGamErr";
  fillParamMapAmps(JPCLSPsiToChiGam, theSuffix, _paramErr.PsiToChiGam);

  std::vector<boost::shared_ptr<const JPCLS> > JPCLSChiToKst1K=psi2SToKpKmPiGamStatesPtr->ChiToKst1KStates();
  theSuffix="K890PiVal";
  fillParamMapAmps(JPCLSChiToKst1K, theSuffix, _paramVal.ChiToK890K);
  theSuffix="K890PiErr";
  fillParamMapAmps(JPCLSChiToKst1K, theSuffix, _paramErr.ChiToK890K);

  std::vector<boost::shared_ptr<const JPCLS> > JPCLSChiToK1400_0_K=psi2SToKpKmPiGamStatesPtr->ChiToKst0KStates();
  theSuffix="K1400_0_PiVal";
  fillParamMapAmps(JPCLSChiToK1400_0_K, theSuffix, _paramVal.ChiToK1400_0_K);
  theSuffix="K1400_0_PiErr";
  fillParamMapAmps(JPCLSChiToK1400_0_K, theSuffix, _paramErr.ChiToK1400_0_K);

  std::vector<boost::shared_ptr<const JPCLS> > JPCLSChiToK1400_1_K=psi2SToKpKmPiGamStatesPtr->ChiToKst1KStates();
  theSuffix="K1400_1_PiVal";
  fillParamMapAmps(JPCLSChiToK1400_1_K, theSuffix, _paramVal.ChiToK1400_1_K);
  theSuffix="K1400_1_PiErr";
  fillParamMapAmps(JPCLSChiToK1400_1_K, theSuffix, _paramErr.ChiToK1400_1_K);

  std::vector<boost::shared_ptr<const JPCLS> > JPCLSChiToK1400_2_K=psi2SToKpKmPiGamStatesPtr->ChiToKst2KStates();
  theSuffix="K1400_2_PiVal";
  fillParamMapAmps(JPCLSChiToK1400_2_K, theSuffix, _paramVal.ChiToK1400_2_K);
  theSuffix="K1400_2_PiErr";
  fillParamMapAmps(JPCLSChiToK1400_2_K, theSuffix, _paramErr.ChiToK1400_2_K);

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToa980pi=psi2SToKpKmPiGamStatesPtr->ChiToa0PiStates();
  theSuffix="a980PiVal";
  fillParamMapAmps(JPCLSChiToa980pi, theSuffix, _paramVal.ChiToa0Pi);
  theSuffix="a980PiErr";
  fillParamMapAmps(JPCLSChiToa980pi, theSuffix, _paramErr.ChiToa0Pi);

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToa2pi=psi2SToKpKmPiGamStatesPtr->ChiToa2PiStates();
  theSuffix="a2PiVal";
  fillParamMapAmps(JPCLSChiToa2pi, theSuffix, _paramVal.ChiToa2Pi);
  theSuffix="a2PiErr";
  fillParamMapAmps(JPCLSChiToa2pi, theSuffix, _paramErr.ChiToa2Pi);

  std::string theKey("K890massVal");
  fillParamMapMass(theKey, _paramVal.BwK890);
  theKey="K890massErr";
  fillParamMapMass(theKey, _paramErr.BwK890);

  theKey="K1400_0massVal";
  fillParamMapMass(theKey, _paramVal.BwK1400_0);
  theKey="K1400_0massErr";
  fillParamMapMass(theKey, _paramErr.BwK1400_0);

  theKey="K1400_1massVal";
  fillParamMapMass(theKey, _paramVal.BwK1400_1);
  theKey="K1400_1massErr";
  fillParamMapMass(theKey, _paramErr.BwK1400_1);

  theKey="K1400_2massVal";
  fillParamMapMass(theKey, _paramVal.BwK1400_2);
  theKey="K1400_2massErr";
  fillParamMapMass(theKey, _paramErr.BwK1400_2);

  theKey="a980massVal";
  fillParamMapMass(theKey, _paramVal.Bwa980);
  theKey="a980massErr";
  fillParamMapMass(theKey, _paramErr.Bwa980);

  theKey="a2massVal";
  fillParamMapMass(theKey, _paramVal.Bwa2);
  theKey="a2massErr";
  fillParamMapMass(theKey, _paramErr.Bwa2);

  theKey="FlatMa980Val";
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

}

void StreamFitParams::fillParamMapAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& toFill){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;

  for ( itJPCLS=theJPCLSs.begin(); itJPCLS!=theJPCLSs.end(); ++itJPCLS){
    //now fill the fitParameterMap
    std::string theKey=(*itJPCLS)->name()+suffix;

    stringPairIter=_stringPairMap.find(theKey);

    if ( stringPairIter != _stringPairMap.end() ){
    std::pair<double, double> theValPair=_stringPairMap[theKey];  
    toFill[*itJPCLS]=theValPair;
    }
    else{
      Warning <<"parameter " << theKey <<" not available in fitVal input file!!!!"  << endmsg;
      std::pair <double,double> theValPair=make_pair(-10000., -10000.);
      toFill[*itJPCLS]=theValPair;      
    }

  }
}

void StreamFitParams::fillParamMapMass(std::string& name, pair<double, double>& toFill){

  StringPairMap::const_iterator stringPairIter;
  stringPairIter=_stringPairMap.find(name);
    if ( stringPairIter != _stringPairMap.end() ){
    std::pair<double, double> theValPair=_stringPairMap[name];  
    toFill=theValPair;
    }
    else{
      Warning <<"parameter " << name <<" not available in fitVal input file!!!!"  << endmsg;
      std::pair <double,double> theValPair=make_pair(-10000., -10000.);
      toFill=theValPair;
    }

}

void StreamFitParams::fillParamFlatte(std::string& name, double& toFill){
  StringPairMap::const_iterator stringPairIter;
  stringPairIter=_stringPairMap.find(name);
  if ( stringPairIter != _stringPairMap.end() ){
    std::pair<double, double> theValPair=_stringPairMap[name];
    toFill=theValPair.first;
  }
  else{
    Warning <<"parameter " << name <<" not available in fitVal input file!!!!"  << endmsg;
    toFill=-10000.;
  }  
}
