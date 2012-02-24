#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp8Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp8Lh::Hyp8Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp7Lh(theEvtList, hypMap )
  ,_K_0_1460ToKstPiHyp8(true)
  ,_K_0_1830ToKstPiHyp8(true)
  ,_K_1_1650Hyp8(true)
  ,_doHyp8(true)
  ,_nFitParams(0)
  ,_evtCounter(0)
  ,_equalParameter(false)
{
  setUp(hypMap); 
}

Hyp8Lh::Hyp8Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp7Lh(theLhPtr->getEventList(), hypMap)
  ,_K_0_1460ToKstPiHyp8(true)
  ,_K_0_1830ToKstPiHyp8(true)
  ,_K_1_1650Hyp8(true)
  ,_doHyp8(true)
  ,_nFitParams(0)
  ,_evtCounter(0)
  ,_equalParameter(false)
{
  setUp(hypMap); 
}

Hyp8Lh::~Hyp8Lh()
{;
}



complex<double> Hyp8Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp7Lh::chi0DecAmps(theParamVal, theData);

  if (!_doHyp8) return result;

  if (_evtCounter==0){
    _equalParameter=equalParams();

    DebugMsg << "equal parameter: "<< _equalParameter << endmsg;

  } 

  if(_equalParameter){
    result+=_currentResultHyp8[_evtCounter];
    _evtCounter++;
    return result;
  }

  complex<double> currentResult(0.,0.);

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;
  double K_0_1430Mass=theParamVal.BwK_0_1430.first;
  double K_0_1430Width=theParamVal.BwK_0_1430.second;

  if (_K_0_1460ToKstPiHyp8){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK892Pi=theParamVal.K_0_1460ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK_0_1430Pi=theParamVal.K_0_1460ToK_0_1430Pi;
    
    double K_0_1460Mass=theParamVal.BwK_0_1460.first;
    double K_0_1460Width=theParamVal.BwK_0_1460.second;
    
    currentResult+=chiToK0KT0K1piKAmp(theData, K_0_1460ToK892Pi, K_0_1460Mass, K_0_1460Width, K892Mass, K892Width);
    currentResult+=chiToK0KT0KpipiKAmp(theData, K_0_1460ToK_0_1430Pi, K_0_1460Mass, K_0_1460Width, K_0_1430Mass, K_0_1430Width); 
  }

  if (_K_0_1830ToKstPiHyp8){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1830ToK892Pi=theParamVal.K_0_1830ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1830ToK_0_1430Pi=theParamVal.K_0_1830ToK_0_1430Pi;
    
    double K_0_1830Mass=theParamVal.BwK_0_1830.first;
    double K_0_1830Width=theParamVal.BwK_0_1830.second;
    
    currentResult+=chiToK0KT0K1piKAmp(theData, K_0_1830ToK892Pi, K_0_1830Mass, K_0_1830Width, K892Mass, K892Width);
    currentResult+=chiToK0KT0KpipiKAmp(theData, K_0_1830ToK_0_1430Pi, K_0_1830Mass, K_0_1830Width, K_0_1430Mass, K_0_1430Width); 
  }

  if (_K_1_1650Hyp8){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1650K=theParamVal.ChiToK_1_1650K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1650ToK892Pi=theParamVal.K_1_1650ToK892Pi;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1650ToK_0_1430Pi=theParamVal.K_1_1650ToK_0_1430Pi;
    
    double K_1_1650Mass=theParamVal.BwK_1_1650.first;
    double K_1_1650Width=theParamVal.BwK_1_1650.second;

    currentResult+=chiToK1ToK1piAmp(theData, ChiToK_1_1650K, K_1_1650ToK892Pi, K_1_1650Mass, K_1_1650Width, K892Mass, K892Width);
    currentResult+=chiToK1ToK0piAmp(theData, K_1_1650ToK_0_1430Pi, K_1_1650ToK_0_1430Pi, K_1_1650Mass, K_1_1650Width, K_0_1430Mass, K_0_1430Width);

  }

  _currentResultHyp8[_evtCounter]=currentResult; 
  _evtCounter++;
  result+=currentResult;
  return result;
}




void Hyp8Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp7Lh::setMnUsrParams(upar, startVal, errVal);

  if (!_doHyp8) return;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){

    _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, (*itAmps));
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, (*itMasses) );  
  }

}



int Hyp8Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }   

  int counter=Hyp7Lh::setFitParamVal(theParamVal, par);

  if (!_doHyp8) return counter;

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

unsigned int  Hyp8Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp7Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp8Lh::print(std::ostream& os) const{
  os << "Hyp8Lh::print\n";
}

void Hyp8Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp7Lh::printCurrentFitResult(theParamVal);

  if (!_doHyp8) return;

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

  if (!_doHyp8) return;

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

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K_0_1460ToKstPiHyp8");

  if (iter !=hypMap.end()){
    _K_0_1460ToKstPiHyp8= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_0_1460ToKstPiHyp8 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "K_0_1460ToKstPiHyp8 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("K_0_1830ToKstPiHyp8");

  if (iter !=hypMap.end()){
    _K_0_1830ToKstPiHyp8= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_0_1830ToKstPiHyp8 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "K_0_1830ToKstPiHyp8 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("K_1_1650Hyp8");

  if (iter !=hypMap.end()){
    _K_1_1650Hyp8= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K_1_1650Hyp8 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "K_1_1650Hyp8 not set!!!" <<endmsg;
    exit(0);
  }
 
  if(!_K_0_1460ToKstPiHyp8 && !_K_0_1830ToKstPiHyp8 && !_K_1_1650Hyp8) _doHyp8=false; 

  if (!_doHyp8) return;

  if (_K_0_1460ToKstPiHyp8){
    _ampVec.push_back(paramEnum2K2PiGam::K_0_1460ToK892Pi);
    _ampVec.push_back(paramEnum2K2PiGam::K_0_1460ToK_0_1430Pi);
    _massVec.push_back(paramEnum2K2PiGam::K_0_1460);
  }

  if (_K_0_1830ToKstPiHyp8){
    _ampVec.push_back(paramEnum2K2PiGam::K_0_1830ToK892Pi);
    _ampVec.push_back(paramEnum2K2PiGam::K_0_1830ToK_0_1430Pi);
    _massVec.push_back(paramEnum2K2PiGam::K_0_1830);
  }

  if (_K_1_1650Hyp8){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToK_1_1650K);
    _ampVec.push_back(paramEnum2K2PiGam::K_1_1650ToK892Pi);
    _ampVec.push_back(paramEnum2K2PiGam::K_1_1650ToK_0_1430Pi);
    _massVec.push_back(paramEnum2K2PiGam::K_1_1650);
  }

  if(!_K0_1430_K0_1430Hyp && !_K0_1430_K0_1430Hyp && !_K1_1270Hyp && !_K0_1430_K892Hyp1) _massVec.push_back(paramEnum2K2PiGam::K_0_1430);

  // fill all other resonances
  if (_K0_1430_K0_1430Hyp || _K0_1430_K0_1430Hyp || _K1_1270Hyp || _K0_1430_K892Hyp1){
  _massVecRemain.push_back(paramEnum2K2PiGam::K_0_1430);
  }

  _massVecRemain.push_back(paramEnum2K2PiGam::K892);

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




void Hyp8Lh::copyCurrentVals(Hyp8Lh* theLh){
  Hyp7Lh::copyCurrentVals(theLh);
  std::map<unsigned int, complex<double> > newResult; 
  std::map<unsigned int, complex<double> >::iterator it;
  for (it= _currentResultHyp8.begin(); it!= _currentResultHyp8.end(); ++it){
    newResult[it->first]=it->second;
  }
  theLh->_currentResultHyp8=newResult;
  
}

bool Hyp8Lh::equalParams(){
  bool result=true;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  if (!compAmpParms( _ampVec )) return false;
  if (!compMassParms(_massVec)) return false;
  if (!compMassParms(_massVecRemain)) return false;

  return result;
}
