#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp9Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp9Lh::Hyp9Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp8Lh(theEvtList, hypMap )
  ,_Pi_2_1670Tof_2_1270PiHyp9(true)
  ,_Pi_2_1670ToK892KHyp9(true)
  ,_doHyp9(true)
  ,_nFitParams(0)
{
  setUp(hypMap); 
}

Hyp9Lh::Hyp9Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp8Lh(theLhPtr->getEventList(), hypMap)
  ,_Pi_2_1670Tof_2_1270PiHyp9(true)
  ,_Pi_2_1670ToK892KHyp9(true)
  ,_doHyp9(true)
  ,_nFitParams(0)
{
  setUp(hypMap); 
}

Hyp9Lh::~Hyp9Lh()
{;
}



complex<double> Hyp9Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp8Lh::chi0DecAmps(theParamVal, theData);

  if (!_doHyp9) return result;
  double Pi_2_1670Mass=theParamVal.BwPi_2_1670.first;
  double Pi_2_1670Width=theParamVal.BwPi_2_1670.second;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi_2_1670Pi=theParamVal.ChiToPi_2_1670Pi; 

  if(_Pi_2_1670Tof_2_1270PiHyp9){
    double f_2_1270Mass=theParamVal.Bwf_2_1270.first;
    double f_2_1270Width=theParamVal.Bwf_2_1270.second;

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_1670Tof_2_1270Pi=theParamVal.Pi_2_1670Tof_2_1270Pi;
    
    result+=chiToPi2Pi0Tof2PiAmp(theData, ChiToPi_2_1670Pi, Pi_2_1670Tof_2_1270Pi, Pi_2_1670Mass, Pi_2_1670Width, f_2_1270Mass, f_2_1270Width);
  
  }

  if(_Pi_2_1670ToK892KHyp9){
    double K892Mass=theParamVal.BwK892.first;
    double K892Width=theParamVal.BwK892.second;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_1670ToK892K=theParamVal.Pi_2_1670ToK892K;

    result+=chiToPi2Pi0ToKstarKAmp(theData, ChiToPi_2_1670Pi, Pi_2_1670ToK892K, Pi_2_1670Mass, Pi_2_1670Width, K892Mass, K892Width);

  }

  return result;
}




void Hyp9Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp8Lh::setMnUsrParams(upar, startVal, errVal);

  if (!_doHyp9) return;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, (*itMasses) );  
  }

}



int Hyp9Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }   

  int counter=Hyp8Lh::setFitParamVal(theParamVal, par);

  if (!_doHyp9) return counter;

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

unsigned int  Hyp9Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp8Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp9Lh::print(std::ostream& os) const{
  os << "Hyp9Lh::print\n";
}

void Hyp9Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp8Lh::printCurrentFitResult(theParamVal);

  if (!_doHyp9) return;

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
    DebugMsg <<"\t mass:" << tmpParam.first <<"\t width:" << tmpParam.second  << endmsg;
  }

}

void Hyp9Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp8Lh::dumpCurrentResult(os, theParamVal, suffix);

  if (!_doHyp9) return;

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

void Hyp9Lh::setUp(const std::map<const std::string, bool>& hypMap){


  if (!_doHyp9) return;

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("Pi_2_1670Tof_2_1270PiHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_1670Tof_2_1270PiHyp9= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _Pi_2_1670Tof_2_1270PiHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_1670Tof_2_1270PiHyp9 not set!!!" <<endmsg;
    exit(0);
  }

  iter= hypMap.find("Pi_2_1670ToK892KHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_1670ToK892KHyp9= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _Pi_2_1670ToK892KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_1670ToK892KHyp9 not set!!!" <<endmsg;
    exit(0);
  }  

  _ampVec.push_back(paramEnum2K2PiGam::ChiToPi_2_1670Pi);
  _massVec.push_back(paramEnum2K2PiGam::Pi_2_1670);

  if (_Pi_2_1670Tof_2_1270PiHyp9){
    _ampVec.push_back(paramEnum2K2PiGam::Pi_2_1670Tof_2_1270Pi);
    _massVec.push_back(paramEnum2K2PiGam::f_2_1270);
  }

  if (_Pi_2_1670ToK892KHyp9){
    _ampVec.push_back(paramEnum2K2PiGam::Pi_2_1670ToK892K);
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
