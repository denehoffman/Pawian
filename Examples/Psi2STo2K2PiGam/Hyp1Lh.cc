#include <getopt.h>
#include <fstream>

#include "Examples/Psi2STo2K2PiGam/Hyp1Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp1Lh::Hyp1Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  AbsPsi2STo2K2PiGamLh(theEvtList)
  ,_K1_1270Hyp(true)
  ,_K1_1400Hyp(false)
  ,_K0_1430_K0_1430Hyp(true)
  ,_K2_1430_K2_1430Hyp(false)
  ,_K0_1430_K2_1430Hyp(false)
  ,_nFitParams(0) 
{
  setUp(hypMap); 
}

Hyp1Lh::Hyp1Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsPsi2STo2K2PiGamLh(theLhPtr->getEventList())
  ,_K1_1270Hyp(true)
  ,_K1_1400Hyp(false)
  ,_K0_1430_K0_1430Hyp(true)
  ,_K2_1430_K2_1430Hyp(false)
  ,_K0_1430_K2_1430Hyp(false)
  ,_nFitParams(0) 
{
  setUp(hypMap);
}

Hyp1Lh::~Hyp1Lh()
{;
}


complex<double> Hyp1Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K892=theParamVal.ChiTo2K892;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_2_1430=theParamVal.ChiTo2K_2_1430;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_0_1430=theParamVal.ChiTo2K_0_1430;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430_K_2_1430=theParamVal.ChiToK_0_1430_K_2_1430;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400ToK892pi=theParamVal.ChiToK1400ToK892pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K1400ToK892Pi=theParamVal.K1400ToK892Pi;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK892Pi=theParamVal.K_1_1270ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK_0_1430Pi=theParamVal.K_1_1270ToK_0_1430Pi;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1710_k=theParamVal.f980_pif1710_k;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1710_pi=theParamVal.f980_kf1710_pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f980=theParamVal.ChiTof980f980;

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;
  double K_2_1430Mass=theParamVal.BwK_2_1430.first;
  double K_2_1430Width=theParamVal.BwK_2_1430.second;
  double K_0_1430Mass=theParamVal.BwK_0_1430.first;
  double K_0_1430Width=theParamVal.BwK_0_1430.second;
  double K1400Mass=theParamVal.BwK_1_1400.first;
  double K1400Width=theParamVal.BwK_1_1400.second;
  double K_1_1270Mass=theParamVal.BwK_1_1270.first;
  double K_1_1270Width=theParamVal.BwK_1_1270.second;
  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  double f1710Mass=theParamVal.Bwf1710.first;
  double f1710Width=theParamVal.Bwf1710.second;



  //Chi_c0 decay to K*(892) K*(892)
  result+=chiTo2K892Amp(theData, ChiTo2K892, K892Mass, K892Width);

  //Chi_c0 decay to K2*(1430) K2*(1430)
  if (_K0_1430_K0_1430Hyp) result+=chiTo2K_0_Amp(theData, ChiTo2K_0_1430, K_0_1430Mass, K_0_1430Width, K_0_1430Mass, K_0_1430Width);
  if (_K2_1430_K2_1430Hyp) result+=chiTo2K_2_Amp(theData, ChiTo2K_2_1430, K_2_1430Mass, K_2_1430Width);
  if (_K0_1430_K2_1430Hyp) result+=chiToK0K2Amp(theData, ChiToK_0_1430_K_2_1430, K_0_1430Mass, K_0_1430Width, K_2_1430Mass, K_2_1430Width);

  //Chi_c0 decay to K1*(1400) -> K1*(892) pi0 -> (K pi0) pi0
  if (_K1_1270Hyp){
    result+=chiToK1ToK1piAmp(theData, K_1_1270ToK892Pi, K_1_1270Mass, K_1_1270Width, K892Mass, K892Width);
    result+=chiToK1ToK0piAmp(theData, K_1_1270ToK_0_1430Pi, K_1_1270Mass, K_1_1270Width, K_0_1430Mass, K_0_1430Width);
  }
  if (_K1_1400Hyp){
    result+=chiToK1ToK1piAmp(theData, K1400ToK892Pi, K1400Mass, K1400Width, K892Mass, K892Width);
  }
  
    //Chi_c0 decay to f0(980) f0(1710) -> (pi0 pi0) (K K) 
  result+=chiTof980_pif0_kAmp(theData, ChiTof980_pif1710_k, f980_Mass, f980_gPiPi,  f980_gKK, f1710Mass, f1710Width);

  //Chi_c0 decay to f0(980) f0(1710) -> (K K) (pi0 pi0) 
  result+=chiTof980_kf0_piAmp(theData, ChiTof980_kf1710_pi, f980_Mass, f980_gKK,  f980_gPiPi, f1710Mass, f1710Width);

  //Chi_c0 decay to f0(980) f0(980) -> (pi0 pi0) (K K) 
  result+=chiTof980f980Amp(theData, ChiTof980f980, f980_Mass, f980_gPiPi,  f980_gKK);

  return result; 

}



void Hyp1Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  _fitParams2K2PiGam.setMnUsrParamsFlatteMass(upar, startVal, errVal, paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace));


  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::ChiGam);
  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K892K892);

  if (_K0_1430_K0_1430Hyp){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1430K_0_1430);
  }
  if (_K2_1430_K2_1430Hyp){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_2_1430K_2_1430);
  }
  if (_K0_1430_K2_1430Hyp){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1430K_2_1430);
  }

  if (_K1_1270Hyp){ 
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_1_1270ToK892Pi);
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_1_1270ToK_0_1430Pi);
  }

  if (_K1_1400Hyp){ 
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_1_1400ToK892Pi);
  }

  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f980_pif1710_k);
  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f980_kf1710_pi);
  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f980f980);

  //K*(892) BW
  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K892);

  if (_K0_1430_K0_1430Hyp || _K0_1430_K2_1430Hyp){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1430);
  }

  if (_K2_1430_K2_1430Hyp || _K0_1430_K2_1430Hyp){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_2_1430);
  }

  //K1(1400) BW
  if (_K1_1270Hyp){ 
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_1_1270);
    if (!_K0_1430_K0_1430Hyp && !_K0_1430_K2_1430Hyp) _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1430); 
  }

  if (_K1_1400Hyp){ 
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_1_1400);
  }

  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::f1710);

//   setMnUsrParamsMass(upar, startVal, errVal, "f980_pi");
  _fitParams2K2PiGam.setMnUsrParamsFlatteMass(upar, startVal, errVal, paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980));

}



int Hyp1Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  int counter=0;
  counter=_fitParams2K2PiGam.setFitParamFlatteMass(theParamVal, par, counter, paramEnum2K2PiGam::name(paramEnum2K2PiGam::phaseSpace));


  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::ChiGam);

  //K*(982) K*(982)    amplitude params
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K892K892);

  //K2*(1430) K2*(1430)    amplitude params
  if (_K0_1430_K0_1430Hyp){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1430K_0_1430);
  }
 if (_K2_1430_K2_1430Hyp){  
   counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_2_1430K_2_1430);
 }

 if (_K0_1430_K2_1430Hyp){  
   counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1430K_2_1430);
 }
  //K1(1400) K    amplitude params
  if (_K1_1270Hyp){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_1_1270ToK892Pi);
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_1_1270ToK_0_1430Pi);
  }

  if (_K1_1400Hyp){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_1_1400ToK892Pi);
  }

  //f980 f1710    amplitude params
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f980_pif1710_k);
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f980_kf1710_pi);

  //f980 f980    amplitude params
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f980f980);

  // K*892 mass 
  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K892);

  if (_K0_1430_K0_1430Hyp || _K0_1430_K2_1430Hyp){
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1430);
  }

  if (_K2_1430_K2_1430Hyp || _K0_1430_K2_1430Hyp){  
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_2_1430);
  }

  //K1(1400) mass
  if (_K1_1270Hyp){
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_1_1270);
    if (!_K0_1430_K0_1430Hyp && !_K0_1430_K2_1430Hyp) counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1430);
  }
  if (_K1_1400Hyp){
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_1_1400);
  }

  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::f1710);

  counter=_fitParams2K2PiGam.setFitParamFlatteMass(theParamVal, par, counter, paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980));


  return counter;
}

unsigned int  Hyp1Lh::nFitParams(){

  return _nFitParams;
}


void Hyp1Lh::print(std::ostream& os) const{
  os << "Hyp1Lh::print\n";
}

void Hyp1Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  DebugMsg<< "phaseSpace: " << theParamVal.phaseSpace << endmsg;


  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);
    
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

    DebugMsg<< "f980:" << endmsg;
    double tmpParamf980Mass=theParamVal.Flatf980;
    double tmpParamf980gPiPi=theParamVal.Flatf980gPiPi;
    double tmpParamf980gKK=theParamVal.Flatf980gKK;
    DebugMsg <<"\t mass:" << tmpParamf980Mass <<"\t gPiPi:" << tmpParamf980gPiPi  <<"\t gKK:" << tmpParamf980gKK << endmsg;



}

void Hyp1Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  std::string tmpStringPs="phaseSpace"+suffix;
  os << tmpStringPs << "\t" << theParamVal.phaseSpace  << "\t" << 0. << std::endl;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);
    
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


  std::string tmpStringFlatMass=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980)+suffix;
  std::string tmpStringFlatPiPi=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gpipi)+suffix;
  std::string tmpStringFlatKK=paramEnum2K2PiGam::name(paramEnum2K2PiGam::f980gKK)+suffix;
  double tmpParamFlatMass = theParamVal.Flatf980;
  double tmpParamFlatPiPi=theParamVal.Flatf980gPiPi;
  double tmpParamFlatKK=theParamVal.Flatf980gKK;
  os << tmpStringFlatMass << "\t" << tmpParamFlatMass <<"\t" << 0.  << std::endl;
  os << tmpStringFlatPiPi << "\t" << tmpParamFlatPiPi <<"\t" << 0.  << std::endl;
  os << tmpStringFlatKK << "\t" << tmpParamFlatKK <<"\t" << 0.  << std::endl;

}


void Hyp1Lh::setUp(const std::map<const std::string, bool>& hypMap){
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K1_1270Hyp");

  if (iter !=hypMap.end()){
    _K1_1270Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1270Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1270Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K1_1400Hyp");
  if (iter !=hypMap.end()){
    _K1_1400Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1400Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1400Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K0_1430_K0_1430Hyp");
  if (iter !=hypMap.end()){
    _K0_1430_K0_1430Hyp= iter->second;
    _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K0_1430_K0_1430Hyp <<endmsg;
  }
  else Alert << "hypothesis K0_1430_K0_1430Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K2_1430_K2_1430Hyp");
  if (iter !=hypMap.end()){
    _K2_1430_K2_1430Hyp= iter->second;
    _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K2_1430_K2_1430Hyp <<endmsg;
  }
  else Alert << "hypothesis K2_1430_K2_1430Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K0_1430_K2_1430Hyp");
  if (iter !=hypMap.end()){
    _K0_1430_K2_1430Hyp= iter->second;
    _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K0_1430_K2_1430Hyp <<endmsg;
  }
  else Alert << "hypothesis K0_1430_K2_1430Hyp not set!!!" <<endmsg;

  _ampVec.push_back(paramEnum2K2PiGam::ChiGam);
  _ampVec.push_back(paramEnum2K2PiGam::K892K892);
  _ampVec.push_back(paramEnum2K2PiGam::f980_pif1710_k);
  _ampVec.push_back(paramEnum2K2PiGam::f980_kf1710_pi);
  _ampVec.push_back(paramEnum2K2PiGam::f980f980);
  if(_K1_1270Hyp){
    _ampVec.push_back(paramEnum2K2PiGam::K_1_1270ToK892Pi);
    _ampVec.push_back(paramEnum2K2PiGam::K_1_1270ToK_0_1430Pi);
  }
  if(_K1_1400Hyp) _ampVec.push_back(paramEnum2K2PiGam::K_1_1400ToK892Pi);
  if(_K2_1430_K2_1430Hyp) _ampVec.push_back(paramEnum2K2PiGam::K_2_1430K_2_1430);
  if(_K0_1430_K0_1430Hyp) _ampVec.push_back(paramEnum2K2PiGam::K_0_1430K_0_1430);
  if(_K0_1430_K2_1430Hyp) _ampVec.push_back(paramEnum2K2PiGam::K_0_1430K_2_1430);

  _massVec.push_back(paramEnum2K2PiGam::K892);
  _massVec.push_back(paramEnum2K2PiGam::f1710);
  if(_K1_1270Hyp)  _massVec.push_back(paramEnum2K2PiGam::K_1_1270);
  if(_K1_1400Hyp)  _massVec.push_back(paramEnum2K2PiGam::K_1_1400);
  if(_K2_1430_K2_1430Hyp || _K0_1430_K2_1430Hyp) _massVec.push_back(paramEnum2K2PiGam::K_2_1430);
  if(_K0_1430_K0_1430Hyp || _K0_1430_K2_1430Hyp || _K1_1270Hyp) _massVec.push_back(paramEnum2K2PiGam::K_0_1430);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  std::vector<unsigned int>::iterator massIt; 
  for (massIt=_massVec.begin(); massIt!=_massVec.end(); ++massIt){
    _nFitParams+=2;
  }

  _nFitParams+=3; //f(980) Flatte parameters
  _nFitParams+=1; //phase space

}
