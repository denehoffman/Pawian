#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp7Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp7Lh::Hyp7Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp6Lh(theEvtList, hypMap )
  ,_KappaHyp(true)
  ,_K1_1680Hyp(true)
  ,_K1_2300Hyp(true)
{
  setUp(hypMap); 
}

Hyp7Lh::Hyp7Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp6Lh(theLhPtr->getEventList(), hypMap)
  ,_KappaHyp(true)
  ,_K1_1680Hyp(true)
  ,_K1_2300Hyp(true)
{
  setUp(hypMap); 
}

Hyp7Lh::~Hyp7Lh()
{;
}



complex<double> Hyp7Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp6Lh::chi0DecAmps(theParamVal, theData);

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;

  //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K

  if(_KappaHyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2Kappa=theParamVal.ChiTo2Kappa;
    double KappaMass=theParamVal.BwKappa.first;
    double KappaWidth=theParamVal.BwKappa.second;
    result+=chiTo2K_0_Amp(theData, ChiTo2Kappa, KappaMass, KappaWidth,  KappaMass, KappaWidth);
  }

  if(_K1_1680Hyp){ 
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K1680=theParamVal.ChiToK892K1680;
    double K_1_1680Mass=theParamVal.BwK_1_1680.first;
    double K_1_1680Width=theParamVal.BwK_1_1680.second;    
    result+=chiToK1K1Amp(theData, ChiToK892K1680, K892Mass, K892Width,  K_1_1680Mass, K_1_1680Width); 
  }

  if(_K1_2300Hyp){ 
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K2300=theParamVal.ChiToK892K2300;
    double K_1_2300Mass=theParamVal.BwK_1_2300.first;
    double K_1_2300Width=theParamVal.BwK_1_2300.second;
    result+=chiToK1K1Amp(theData, ChiToK892K2300, K892Mass, K892Width,  K_1_2300Mass, K_1_2300Width); 
  }

  return result;
}




void Hyp7Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp6Lh::setMnUsrParams(upar, startVal, errVal);

  if(_KappaHyp){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::KappaKappa);
    if (!_KappaK_0_1950Hyp6) _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::Kappa);
  }

  if(_K1_1680Hyp){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K892K_1_1680);
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_1_1680);
  }

  if(_K1_2300Hyp){
    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::K892K_1_2300);
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::K_1_2300);
  }

}



int Hyp7Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }   

  int counter=Hyp6Lh::setFitParamVal(theParamVal, par);
  if(_KappaHyp){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::KappaKappa);
    if (!_KappaK_0_1950Hyp6) counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::Kappa);  
  }
  if(_K1_1680Hyp){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K892K_1_1680);
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_1_1680);  
  }
  if(_K1_2300Hyp){
    counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::K892K_1_2300);
    counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::K_1_2300);  
  }


  return counter;
}

unsigned int  Hyp7Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp6Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp7Lh::print(std::ostream& os) const{
  os << "Hyp7Lh::print\n";
}

void Hyp7Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp6Lh::printCurrentFitResult(theParamVal);

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

void Hyp7Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp6Lh::dumpCurrentResult(os, theParamVal, suffix);

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

void Hyp7Lh::setUp(const std::map<const std::string, bool>& hypMap){

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("KappaHyp7");
  
  if (iter !=hypMap.end()){
    _KappaHyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _KappaHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis KappaHyp7 not set!!!" <<endmsg;
  
  iter= hypMap.find("K1_1680Hyp7");
  if (iter !=hypMap.end()){
    _K1_1680Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1680Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1680Hyp7 not set!!!" <<endmsg;
  
  iter= hypMap.find("K1_2300Hyp7");
  
  if (iter !=hypMap.end()){
    _K1_2300Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_2300Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_2300Hyp7 not set!!!" <<endmsg;


  if(_KappaHyp){
    _ampVec.push_back(paramEnum2K2PiGam::KappaKappa);

    if(!_KappaK_0_1950Hyp6) _massVec.push_back(paramEnum2K2PiGam::Kappa);
  }

  if(_K1_1680Hyp){
    _ampVec.push_back(paramEnum2K2PiGam::K892K_1_1680);
    _massVec.push_back(paramEnum2K2PiGam::K_1_1680);
  }

  if(_K1_2300Hyp){
    _ampVec.push_back(paramEnum2K2PiGam::K892K_1_2300);
    _massVec.push_back(paramEnum2K2PiGam::K_1_2300);
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
