#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp6Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp6Lh::Hyp6Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap, bool cacheAmps ) :
  Hyp5Lh(theEvtList, hypMap, cacheAmps )
  ,_doHyp6(true)
  ,_K_0_1430K_0_1950Hyp6(true)
  ,_KappaK_0_1430Hyp6(true)
  ,_KappaK_0_1950Hyp6(true) 
  ,_nFitParams(0)
  ,_evtCounter(0)
  ,_equalParameter(false) 
{
  setUp(hypMap); 
}

Hyp6Lh::Hyp6Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap, bool cacheAmps ) :
  Hyp5Lh(theLhPtr->getEventList(), hypMap, cacheAmps)
  ,_doHyp6(true)
  ,_K_0_1430K_0_1950Hyp6(true)
  ,_KappaK_0_1430Hyp6(true)
  ,_KappaK_0_1950Hyp6(true)  
  ,_nFitParams(0)
  ,_evtCounter(0)
  ,_equalParameter(false) 
{
  setUp(hypMap); 
}

Hyp6Lh::~Hyp6Lh()
{;
}

bool  Hyp6Lh::equalChic0DecParams(){
  bool result=false; 
  bool equalRemainHyps=Hyp5Lh::equalChic0DecParams();
  if(!_doHyp6) return equalRemainHyps;
  _equalParameter=equalParams();
  DebugMsg << "equal parameter: "<< _equalParameter << endmsg;

  if(_equalParameter && equalRemainHyps) result=true;
  return result;
}

complex<double> Hyp6Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp5Lh::chi0DecAmps(theParamVal, theData);

  if(!_doHyp6) return result;

  if(_equalParameter){
    result+=_currentResultHyp6[_evtCounter];
    _evtCounter++;
    return result;
  }

  complex<double> currentResult(0.,0.);


  double K_0_1950Mass=theParamVal.BwK_0_1950.first;
  double K_0_1950Width=theParamVal.BwK_0_1950.second;

  if (_K_0_1430K_0_1950Hyp6){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430K_0_1950=theParamVal.ChiToK_0_1430K_0_1950;
    double K_0_1430Mass=theParamVal.BwK_0_1430.first;
    double K_0_1430Width=theParamVal.BwK_0_1430.second;
    
    //Chi_c0 decay to K_0_1950 K_0_1430 
    currentResult+=chiTo2K_0_Amp(theData, ChiToK_0_1430K_0_1950, K_0_1430Mass, K_0_1430Width, K_0_1950Mass, K_0_1950Width); 
  }

  if (_KappaK_0_1430Hyp6){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappaK_0_1430=theParamVal.ChiToKappaK_0_1430;
    double K_0_1430Mass=theParamVal.BwK_0_1430.first;
    double K_0_1430Width=theParamVal.BwK_0_1430.second;
    double KappaMass=theParamVal.BwKappa.first;
    double KappaWidth=theParamVal.BwKappa.second;
    currentResult+=chiTo2K_0_Amp(theData, ChiToKappaK_0_1430, KappaMass, KappaWidth, K_0_1430Mass, K_0_1430Width); 
  }

  if (_KappaK_0_1950Hyp6){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappaK_0_1950=theParamVal.ChiToKappaK_0_1950;
    double KappaMass=theParamVal.BwKappa.first;
    double KappaWidth=theParamVal.BwKappa.second;
    currentResult+=chiTo2K_0_Amp(theData, ChiToKappaK_0_1950, KappaMass, KappaWidth, K_0_1950Mass, K_0_1950Width); 
  }

  _currentResultHyp6[_evtCounter]=currentResult; 
  _evtCounter++;
  result+=currentResult;
  return result;
}




void Hyp6Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal, param2K2PiGam& errVal){

  Hyp5Lh::setMnUsrParams(upar, startVal, errVal);

  if(!_doHyp6) return;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, (*itMasses) );  
  }

}



int Hyp6Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }  
  
  int counter=Hyp5Lh::setFitParamVal(theParamVal, par);

  if(!_doHyp6) return counter;

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

unsigned int  Hyp6Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp5Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp6Lh::print(std::ostream& os) const{
  os << "Hyp6Lh::print\n";
}

void Hyp6Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp5Lh::printCurrentFitResult(theParamVal);
  if(!_doHyp6) return;
  
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

void Hyp6Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp5Lh::dumpCurrentResult(os, theParamVal, suffix);

  if(!_doHyp6) return;

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

void Hyp6Lh::setUp(const std::map<const std::string, bool>& hypMap){


  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K_0_1430K_0_1950Hyp6");

  if (iter !=hypMap.end()){
    _K_0_1430K_0_1950Hyp6= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_0_1430K_0_1950Hyp6 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K_0_1430K_0_1950Hyp6 not set!!!" <<endmsg; 

  iter= hypMap.find("KappaK_0_1430Hyp6");

  if (iter !=hypMap.end()){
    _KappaK_0_1430Hyp6= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _KappaK_0_1430Hyp6 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis KappaK_0_1430Hyp6 not set!!!" <<endmsg; 


  iter= hypMap.find("KappaK_0_1950Hyp6");

  if (iter !=hypMap.end()){
    _KappaK_0_1950Hyp6= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _KappaK_0_1950Hyp6 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis KappaK_0_1950Hyp6 not set!!!" <<endmsg; 

  if(!_K_0_1430K_0_1950Hyp6 && !_KappaK_0_1950Hyp6 && !_KappaK_0_1430Hyp6){
    _doHyp6=false;
    return;
  }


  if (_KappaK_0_1430Hyp6){
  _ampVec.push_back(paramEnum2K2PiGam::ChiToKappaK_0_1430);
  }

  if(_K_0_1430K_0_1950Hyp6){
  _ampVec.push_back(paramEnum2K2PiGam::K_0_1430K_0_1950);

  _massVec.push_back(paramEnum2K2PiGam::K_0_1950);
  }

  if (_K_0_1430K_0_1950Hyp6 || _KappaK_0_1430Hyp6){
    if(!_K0_1430_K0_1430Hyp && !_K0_1430_K2_1430Hyp && !_K1_1270Hyp && !_K0_1430_K892Hyp1) _massVec.push_back(paramEnum2K2PiGam::K_0_1430);
  }


  if(_KappaK_0_1950Hyp6){
  _ampVec.push_back(paramEnum2K2PiGam::KappaK_0_1950);

  _massVec.push_back(paramEnum2K2PiGam::Kappa);
  if(!_K_0_1430K_0_1950Hyp6) _massVec.push_back(paramEnum2K2PiGam::K_0_1950);
  }

  // fill all other resonances
  if (_K_0_1430K_0_1950Hyp6 || _KappaK_0_1430Hyp6){
    if (_K0_1430_K0_1430Hyp || _K0_1430_K0_1430Hyp || _K1_1270Hyp || _K0_1430_K892Hyp1){
      _massVecRemain.push_back(paramEnum2K2PiGam::K_0_1430);
    }
  }
  if(_KappaK_0_1430Hyp6 || _KappaK_0_1950Hyp6){
    _massVecRemain.push_back(paramEnum2K2PiGam::Kappa);
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

void Hyp6Lh::copyCurrentVals(Hyp6Lh* theLh){
  Hyp5Lh::copyCurrentVals(theLh);
  std::map<unsigned int, complex<double> > newResult; 
  std::map<unsigned int, complex<double> >::iterator it;
  for (it= _currentResultHyp6.begin(); it!= _currentResultHyp6.end(); ++it){
    newResult[it->first]=it->second;
  }
  theLh->_currentResultHyp6=newResult;
  
}

bool Hyp6Lh::equalParams(){
  bool result=true;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  if (!compAmpParms( _ampVec )) return false;
  if (!compMassParms(_massVec)) return false;
  if (!compMassParms(_massVecRemain)) return false;

  return result;
}
