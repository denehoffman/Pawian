#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp4Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp4Lh::Hyp4Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp3Lh(theEvtList, hypMap )
  ,_f980f1370Hyp4(true)
  ,_f1710f1370Hyp4(true)
  ,_disableHyp4(false)
  ,_nFitParams(0)
{
  setUp(hypMap);
}

Hyp4Lh::Hyp4Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap  ) :
  Hyp3Lh(theLhPtr->getEventList(), hypMap)
  ,_f980f1370Hyp4(true)
  ,_f1710f1370Hyp4(true)
  ,_disableHyp4(false)
  ,_nFitParams(0)
{
  setUp(hypMap);
}

Hyp4Lh::~Hyp4Lh()
{;
}



complex<double> Hyp4Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp3Lh::chi0DecAmps(theParamVal, theData);

  if(_disableHyp4) return result;
  double f1370Mass=theParamVal.Bwf1370.first;
  double f1370Width=theParamVal.Bwf1370.second;

  if (_f980f1370Hyp4){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1370_k=theParamVal.f980_pif1370_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1370_pi=theParamVal.f980_kf1370_pi;
    double f980_Mass=theParamVal.Flatf980;
    double f980_gPiPi=theParamVal.Flatf980gPiPi;
    double f980_gKK=theParamVal.Flatf980gKK;
    
    //Chi_c0 decay to f0(980) f0(1370) -> (pi0 pi0) (K K) 
    result+=chiTof980_pif0_kAmp(theData, ChiTof980_pif1370_k, f980_Mass, f980_gPiPi,  f980_gKK, f1370Mass, f1370Width);
    result+=chiTof980_kf0_piAmp(theData, ChiTof980_kf1370_pi, f980_Mass, f980_gKK,  f980_gPiPi, f1370Mass, f1370Width);
  }

  if (_f1710f1370Hyp4){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof1710_pif1370_k=theParamVal.f1710_pif1370_k;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof1710_kf1370_pi=theParamVal.f1710_kf1370_pi;
    double f1710Mass=theParamVal.Bwf1710.first;
    double f1710Width=theParamVal.Bwf1710.second;

    result+=chiTof0_pif0_kAmp(theData, ChiTof1710_pif1370_k, f1710Mass, f1710Width,  f1370Mass, f1370Width);
    result+=chiTof0_pif0_kAmp(theData, ChiTof1710_kf1370_pi, f1370Mass, f1370Width,  f1710Mass, f1710Width);
  }

  return result;
}


void Hyp4Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal, param2K2PiGam& errVal){

  Hyp3Lh::setMnUsrParams(upar, startVal, errVal);
  if (_disableHyp4) return;

  if (_f980f1370Hyp4){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f980_pif1370_k);
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f980_kf1370_pi);
  }

  if (_f1710f1370Hyp4){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f1710_pif1370_k);
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f1710_kf1370_pi);
  }

  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::f1370);  
}



int Hyp4Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }  

  int counter=Hyp3Lh::setFitParamVal(theParamVal, par);

  if (_disableHyp4) return counter;
 
 //sigma f1710    amplitude params
  if (_f980f1370Hyp4){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f980_pif1370_k);
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f980_kf1370_pi);
  }

  if (_f1710f1370Hyp4){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f1710_pif1370_k);
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f1710_kf1370_pi);
  }

  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::f1370);  

  return counter;
}

unsigned int  Hyp4Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp3Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp4Lh::print(std::ostream& os) const{
  os << "Hyp4Lh::print\n";
}

void Hyp4Lh::printCurrentFitResult(param2K2PiGam& theParamVal){
  Hyp3Lh::printCurrentFitResult(theParamVal);

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

void Hyp4Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }
  Hyp3Lh::dumpCurrentResult(os, theParamVal, suffix);

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

void Hyp4Lh::setUp(const std::map<const std::string, bool>& hypMap){

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("f980f1370Hyp4");

  if (iter !=hypMap.end()){
    _f980f1370Hyp4= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f980f1370Hyp4 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "hypothesis f980f1370Hyp4 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("f1710f1370Hyp4");

  if (iter !=hypMap.end()){
    _f1710f1370Hyp4= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _f1710f1370Hyp4 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "hypothesis f1710f1370Hyp4 not set!!!" <<endmsg;
    exit(0);
  }

  if(!_f980f1370Hyp4 && !_f1710f1370Hyp4) _disableHyp4=true; 

  if (_disableHyp4) return;

  if (_f980f1370Hyp4){
    _ampVec.push_back(paramEnum2K2PiGam::f980_pif1370_k);
    _ampVec.push_back(paramEnum2K2PiGam::f980_kf1370_pi);
  }
  if (_f1710f1370Hyp4){
    _ampVec.push_back(paramEnum2K2PiGam::f1710_pif1370_k);
    _ampVec.push_back(paramEnum2K2PiGam::f1710_kf1370_pi);
  }

  _massVec.push_back(paramEnum2K2PiGam::f1370);

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
