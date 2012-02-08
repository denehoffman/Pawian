#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp5Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp5Lh::Hyp5Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp4Lh(theEvtList, hypMap )
  , _doHyp5(true)
  , _K_0_2400KHyp5(true)
  ,_K_0_2400KTof_0_1710Hyp5(true)
  ,_K_1_2400KHyp5(true)
  ,_ChiToK_0_1430KPiHyp5(true)
  ,_ChiToK892KPiHyp5(true)
  , _nFitParams(0)
{
  setUp(hypMap); 
}

Hyp5Lh::Hyp5Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp4Lh(theLhPtr->getEventList(), hypMap)
  , _doHyp5(true)
  , _K_0_2400KHyp5(true)
  ,_K_0_2400KTof_0_1710Hyp5(true)
  ,_K_1_2400KHyp5(true)
  ,_ChiToK_0_1430KPiHyp5(true)
  ,_ChiToK892KPiHyp5(true)
  , _nFitParams(0)
{
  setUp(hypMap); 
}

Hyp5Lh::~Hyp5Lh()
{;
}



complex<double> Hyp5Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp4Lh::chi0DecAmps(theParamVal, theData);
  if (!_doHyp5) return result;

  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  if (_K_0_2400KHyp5 || _K_0_2400KTof_0_1710Hyp5){

    double K_0_2400Mass=theParamVal.BwK_0_2400.first;
    double K_0_2400Width=theParamVal.BwK_0_2400.second;

    if (_K_0_2400KHyp5){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_2400ToKf980=theParamVal.ChiToK_0_2400ToKf980; 
    //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K
    result+=chiToK_0_KToKf980KAmp(theData, ChiToK_0_2400ToKf980, f980_Mass, f980_gPiPi,  f980_gKK, K_0_2400Mass, K_0_2400Width);
    }
    if (_K_0_2400KTof_0_1710Hyp5){
    double f1710_Mass=theParamVal.Bwf1710.first;
    double f1710_Width=theParamVal.Bwf1710.second;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_2400ToKf1710=theParamVal.ChiToK_0_2400ToKf_0_1710; 
    //Chi_c0 decay to K_0_2400 K  -> (K f0(1710) ) K -> (K pi0 pi0 ) K
    result+=chiToK_0_KToKf0KAmp(theData, ChiToK_0_2400ToKf1710, f1710_Mass, f1710_Width, K_0_2400Mass, K_0_2400Width);
    }
 
  }

  if (_K_1_2400KHyp5){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_2400K=theParamVal.ChiToK_1_2400K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_2400Tof980K=theParamVal.K_1_2400Tof980K;
    double K_1_2400Mass=theParamVal.BwK_1_2400.first;
    double K_1_2400Width=theParamVal.BwK_1_2400.second;
    result+=chiToK1Tof980_piKAmp(theData, ChiToK_1_2400K, K_1_2400Tof980K, K_1_2400Mass, K_1_2400Width, f980_Mass, f980_gKK, f980_gPiPi); 
  }

  if(_ChiToK_0_1430KPiHyp5){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430KPi=theParamVal.ChiToK_0_1430KPi;
    double K_0_1430Mass=theParamVal.BwK_0_1430.first;
    double K_0_1430Width=theParamVal.BwK_0_1430.second;
    result+=chiToK0KpiAmp(theData, ChiToK_0_1430KPi, K_0_1430Mass, K_0_1430Width); 

  }

  if(_ChiToK892KPiHyp5){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892KPi=theParamVal.ChiToK892KPi;
    double K892Mass=theParamVal.BwK892.first;
    double K892Width=theParamVal.BwK892.second;
    result+=chiToKst1KpiAmp(theData, ChiToK892KPi, K892Mass, K892Width);
  }
  return result;
}




void Hyp5Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp4Lh::setMnUsrParams(upar, startVal, errVal);

  if (!_doHyp5) return;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, (*itMasses) );  
  }

}



int Hyp5Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){
  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }   

  int counter=Hyp4Lh::setFitParamVal(theParamVal, par);

  if (!_doHyp5) return counter;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, (*itAmps));
  } 

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, (*itMasses) ); 
  }


  return counter;
}

unsigned int  Hyp5Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp4Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp5Lh::print(std::ostream& os) const{
  os << "Hyp5Lh::print\n";
}

void Hyp5Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp4Lh::printCurrentFitResult(theParamVal);

  if (!_doHyp5) return;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;    
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< paramEnum2K2PiGam::name(*itAmps) << endmsg;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParams2K2PiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }  
  }
  
  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    DebugMsg<< paramEnum2K2PiGam::name(*itMasses) << endmsg;
    std::pair<double, double> tmpParam=_fitParams2K2PiGam.massPair(theParamVal, *itMasses);
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

}

void Hyp5Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp4Lh::dumpCurrentResult(os, theParamVal, suffix);

  if (!_doHyp5) return;
  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);

    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;    
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParams2K2PiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];

      std::string tmpStringDec=(*itJPCLS)->name()+paramEnum2K2PiGam::name(*itAmps)+suffix;
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }  
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    std::string tmpStringMass=paramEnum2K2PiGam::name(*itMasses)+suffix;

    std::pair<double, double> tmpParam=_fitParams2K2PiGam.massPair(theParamVal, *itMasses);
    os << tmpStringMass << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

}

void Hyp5Lh::setUp(const std::map<const std::string, bool>& hypMap){

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("doHyp5");

  if (iter !=hypMap.end()){
    _doHyp5= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _doHyp5 <<endmsg;
  }
  else { Alert << "hypothesis disableHyp5 not set!!!" <<endmsg;
    exit(0);
  }


  if (!_doHyp5) return;


  iter= hypMap.find("K_0_2400KHyp5");

  if (iter !=hypMap.end()){
    _K_0_2400KHyp5= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_0_2400KHyp5 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K_0_2400KHyp5 not set!!!" <<endmsg;


  iter= hypMap.find("K_1_2400KHyp5");

  if (iter !=hypMap.end()){
    _K_1_2400KHyp5= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_1_2400KHyp5 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K_1_2400KHyp5 not set!!!" <<endmsg;

  iter= hypMap.find("K_0_2400KTof_0_1710Hyp5");

  if (iter !=hypMap.end()){
    _K_0_2400KTof_0_1710Hyp5= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_0_2400KTof_0_1710Hyp5 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K_0_2400KTof_0_1710Hyp5 not set!!!" <<endmsg;

  iter= hypMap.find("ChiToK_0_1430KPiHyp5");

  if (iter !=hypMap.end()){
    _ChiToK_0_1430KPiHyp5= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _ChiToK_0_1430KPiHyp5 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis ChiToK_0_1430KPiHyp5 not set!!!" <<endmsg;

  iter= hypMap.find("ChiToK892KPiHyp5");
  if (iter !=hypMap.end()){
    _ChiToK892KPiHyp5= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _ChiToK892KPiHyp5 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis ChiToK892KPiHyp5 not set!!!" <<endmsg;

  if(_K_0_2400KHyp5 || _K_0_2400KTof_0_1710Hyp5){

    _massVec.push_back(paramEnum2K2PiGam::K_0_2400);

    if (_K_0_2400KHyp5){
      _ampVec.push_back(paramEnum2K2PiGam::K_0_2400KToKf980);
    }

    if (_K_0_2400KTof_0_1710Hyp5){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToK_0_2400ToKf_0_1710);
    }

  }

  if(_K_1_2400KHyp5){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToK_1_2400K);
    _ampVec.push_back(paramEnum2K2PiGam::K_1_2400Tof980K);
    _massVec.push_back(paramEnum2K2PiGam::K_1_2400);
  }

  if(_ChiToK_0_1430KPiHyp5){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToK_0_1430KPi);
  }

  if (_ChiToK892KPiHyp5){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToK892KPi);
  }

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  std::vector<unsigned int>::iterator massIt; 
  for (massIt=_massVec.begin(); massIt!=_massVec.end(); ++massIt){
    _nFitParams+=2;
  }
}
