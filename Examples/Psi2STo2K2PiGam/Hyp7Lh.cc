#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp7Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp7Lh::Hyp7Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps ) :
  Hyp6Lh(theEvtList, hypMap, theStatesPtr, cacheAmps )
  ,_KappaHyp(true)
  ,_K1_1680Hyp(true)
  ,_K1_1680K1_1680Hyp7(true)
  ,_K1_1680K0_1430Hyp7(true)
  ,_K1_2300Hyp(true)
  ,_nFitParams(0)
  ,_equalParameter(false)
  ,_equalK1_1680Params(false)
  ,_equalK1_2300Params(false)
{
  setUp(hypMap); 
}

Hyp7Lh::Hyp7Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps ) :
  Hyp6Lh(theLhPtr->getEventList(), hypMap, theStatesPtr, cacheAmps)
  ,_KappaHyp(true)
  ,_K1_1680Hyp(true)
  ,_K1_1680K1_1680Hyp7(true)
  ,_K1_1680K0_1430Hyp7(true)
  ,_K1_2300Hyp(true)
  ,_nFitParams(0)
  ,_equalParameter(false)
  ,_equalK1_1680Params(false)
  ,_equalK1_2300Params(false)
{
  setUp(hypMap); 
}

Hyp7Lh::~Hyp7Lh()
{;
}

bool  Hyp7Lh::equalChic0DecParams(){
  bool result=false; 
  bool equalRemainHyps=Hyp6Lh::equalChic0DecParams();
  _equalParameter=equalParams();
  if (compAmpParms(_ampVecK1_1680) && compMassParms(_massVecK1_1680)) _equalK1_1680Params=true;
  if (compAmpParms(_ampVecK1_2300) && compMassParms(_massVecK1_2300)) _equalK1_2300Params=true;
  DebugMsg << "equal parameter: "<< _equalParameter << endmsg;
  if(_equalParameter && equalRemainHyps) result=true;
  return result;
}

complex<double> Hyp7Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp6Lh::chi0DecAmps(theParamVal, theData);

  if(_equalParameter){
    result+=_currentResultHyp7[_evtCounter];
    return result;
  }

  complex<double> currentResult(0.,0.);



  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;

  if(_KappaHyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2Kappa=theParamVal.ChiTo2Kappa;
    double KappaMass=theParamVal.BwKappa.first;
    double KappaWidth=theParamVal.BwKappa.second;
    currentResult+=chiTo2K_0_Amp(theData, ChiTo2Kappa, KappaMass, KappaWidth,  KappaMass, KappaWidth);
  }

  if(_K1_1680Hyp || _K1_1680K1_1680Hyp7 || _K1_1680K0_1430Hyp7){ 

    double K_1_1680Mass=theParamVal.BwK_1_1680.first;
    double K_1_1680Width=theParamVal.BwK_1_1680.second; 

    if (_K1_1680Hyp){
      if(_equalK1_1680Params) currentResult+=_currentResultK1_1680Hyp7[_evtCounter];
      else{
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K1680=theParamVal.ChiToK892K1680;
      complex<double> tmpK1_1680Amp=chiToK1K1Amp(theData, ChiToK892K1680, K892Mass, K892Width,  K_1_1680Mass, K_1_1680Width);
      _currentResultK1_1680Hyp7[_evtCounter]=tmpK1_1680Amp;
      currentResult+=tmpK1_1680Amp;
      }
    } 

    if (_K1_1680K1_1680Hyp7){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1680K1680=theParamVal.ChiToK1680K1680;
      currentResult+=chiToK1K1Amp(theData, ChiToK1680K1680, K_1_1680Mass, K_1_1680Width,  K_1_1680Mass, K_1_1680Width);
    } 

    if(_K1_1680K0_1430Hyp7){

      double K_0_1430Mass=theParamVal.BwK_0_1430.first;
      double K_0_1430Width=theParamVal.BwK_0_1430.second;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1680K_0_1430=theParamVal.ChiToK1680K_0_1430;
      currentResult+=chiToK0K1Amp(theData, ChiToK1680K_0_1430, K_0_1430Mass, K_0_1430Width,  K_1_1680Mass, K_1_1680Width);
    }

  }

  if(_K1_2300Hyp){
    if(_equalK1_2300Params) currentResult+=_currentResultK1_2300Hyp7[_evtCounter];
    else{
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K2300=theParamVal.ChiToK892K2300;
      double K_1_2300Mass=theParamVal.BwK_1_2300.first;
      double K_1_2300Width=theParamVal.BwK_1_2300.second;
      complex<double> tmpK1_2300Amp=chiToK1K1Amp(theData, ChiToK892K2300, K892Mass, K892Width,  K_1_2300Mass, K_1_2300Width);
      _currentResultK1_2300Hyp7[_evtCounter]=tmpK1_2300Amp;
      currentResult+=tmpK1_2300Amp;
    } 
  }

  if(_cacheAmps) _currentResultHyp7[_evtCounter]=currentResult; 
  result+=currentResult;
  return result;
}




void Hyp7Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp6Lh::setMnUsrParams(upar, startVal, errVal);

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, (*itMasses) );  
  }

}



int Hyp7Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }   

  int counter=Hyp6Lh::setFitParamVal(theParamVal, par);

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

  std::string theKey="KappaHyp7";
  setHyps( hypMap, _KappaHyp, theKey);

  theKey="K1_1680Hyp7";
  setHyps( hypMap, _K1_1680Hyp, theKey);

  theKey="K1_1680K1_1680Hyp7";
  setHyps( hypMap, _K1_1680K1_1680Hyp7, theKey);  

  theKey="K1_1680K0_1430Hyp7";
  setHyps( hypMap, _K1_1680K0_1430Hyp7, theKey); 

  theKey="K1_2300Hyp7";
  setHyps( hypMap, _K1_2300Hyp, theKey);


  if(_KappaHyp){
    _ampVec.push_back(paramEnum2K2PiGam::KappaKappa);

    if(!_KappaK_0_1950Hyp6) _massVec.push_back(paramEnum2K2PiGam::Kappa);
  }

  if(_K1_1680Hyp || _K1_1680K1_1680Hyp7 || _K1_1680K0_1430Hyp7){

    _massVec.push_back(paramEnum2K2PiGam::K_1_1680);
    if(_K1_1680Hyp){
      _ampVec.push_back(paramEnum2K2PiGam::K892K_1_1680);
      _ampVecK1_1680.push_back(paramEnum2K2PiGam::K892K_1_1680);
      _massVecK1_1680.push_back(paramEnum2K2PiGam::K892);
      _massVecK1_1680.push_back(paramEnum2K2PiGam::K_1_1680);
    }
    if(_K1_1680K1_1680Hyp7){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToK1680K1680);
    }

    if(_K1_1680K0_1430Hyp7){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToK1680K_0_1430);
    }

  }

  if(_K1_2300Hyp){
    _ampVec.push_back(paramEnum2K2PiGam::K892K_1_2300);
    _massVec.push_back(paramEnum2K2PiGam::K_1_2300);
    _ampVecK1_2300.push_back(paramEnum2K2PiGam::K892K_1_2300);
    _massVecK1_2300.push_back(paramEnum2K2PiGam::K_1_2300);
    _massVecK1_2300.push_back(paramEnum2K2PiGam::K892);
  }

  // fill all other resonances
  if(_KappaK_0_1950Hyp6){
    _massVecRemain.push_back(paramEnum2K2PiGam::Kappa);
  }

  if(_K1_1680Hyp || _K1_2300Hyp){
    _massVecRemain.push_back(paramEnum2K2PiGam::K892);
  }

  if(_K1_1680K0_1430Hyp7){
    _massVecRemain.push_back(paramEnum2K2PiGam::K_0_1430);
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



void Hyp7Lh::copyCurrentVals(Hyp7Lh* theLh){
  Hyp6Lh::copyCurrentVals(theLh);
  if(_cacheAmps){
    theLh->_currentResultHyp7=_currentResultHyp7;
    theLh->_currentResultK1_1680Hyp7=_currentResultK1_1680Hyp7;
    theLh->_currentResultK1_2300Hyp7=_currentResultK1_2300Hyp7;
  }  
}



bool Hyp7Lh::equalParams(){

  bool result=true;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  if (!compAmpParms( _ampVec )) return false;
  if (!compMassParms(_massVec)) return false;
  if (!compMassParms(_massVecRemain)) return false;

  return result;

}
