#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp8Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp8Lh::Hyp8Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp7Lh(theEvtList, hypMap )
{
  setUp(hypMap); 
}

Hyp8Lh::Hyp8Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp7Lh(theLhPtr->getEventList(), hypMap)
{
  setUp(hypMap); 
}

Hyp8Lh::~Hyp8Lh()
{;
}



complex<double> Hyp8Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp7Lh::chi0DecAmps(theParamVal, theData);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK892Pi=theParamVal.K_0_1460ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK_0_1430Pi=theParamVal.K_0_1460ToK_0_1430Pi;

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;
  double K_0_1430Mass=theParamVal.BwK_0_1430.first;
  double K_0_1430Width=theParamVal.BwK_0_1430.second;
  double K_0_1460Mass=theParamVal.BwK_0_1460.first;
  double K_0_1460Width=theParamVal.BwK_0_1460.second;

  //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K


  result+=chiToK0KT0K1piKAmp(theData, K_0_1460ToK892Pi, K_0_1460Mass, K_0_1460Width, K892Mass, K892Width);
  result+=chiToK0KT0KpipiKAmp(theData, K_0_1460ToK_0_1430Pi, K_0_1460Mass, K_0_1460Width, K_0_1430Mass, K_0_1430Width); 


  return result;
}




void Hyp8Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp7Lh::setMnUsrParams(upar, startVal, errVal);

  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1460ToK892Pi);
  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1460ToK_0_1430Pi);
  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1460);
  if(!_K0_1430_K0_1430Hyp && !_K0_1430_K0_1430Hyp && !_K1_1270Hyp) _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_0_1430);
}



int Hyp8Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }   

  int counter=Hyp7Lh::setFitParamVal(theParamVal, par);

  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1460ToK892Pi);
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1460ToK_0_1430Pi);
  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1460);  
  if(!_K0_1430_K0_1430Hyp && !_K0_1430_K0_1430Hyp && !_K1_1270Hyp) counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_0_1430);  
  return counter;
}

unsigned int  Hyp8Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp7Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp8Lh::print(std::ostream& os) const{
  os << "Hyp8Lh::print\n";
}

void Hyp8Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp7Lh::printCurrentFitResult(theParamVal);

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

void Hyp8Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp7Lh::dumpCurrentResult(os, theParamVal, suffix);

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

void Hyp8Lh::setUp(const std::map<const std::string, bool>& hypMap){
  _ampVec.push_back(paramEnum2K2PiGam::K_0_1460ToK892Pi);
  _ampVec.push_back(paramEnum2K2PiGam::K_0_1460ToK_0_1430Pi);

  _massVec.push_back(paramEnum2K2PiGam::K_0_1460);
  if(!_K0_1430_K0_1430Hyp && !_K0_1430_K0_1430Hyp && !_K1_1270Hyp) _massVec.push_back(paramEnum2K2PiGam::K_0_1430);


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
