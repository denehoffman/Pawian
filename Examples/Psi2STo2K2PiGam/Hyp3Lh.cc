#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp3Lh::Hyp3Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps ) :
  Hyp2Lh(theEvtList, hypMap, theStatesPtr, cacheAmps)
  ,_sigmaf980Hyp(false)
  ,_sigmaf1710Hyp(false)
  ,_sigmaf2200Hyp(false)  
  ,_doHyp3(true)
  ,_nFitParams(0)
  ,_f980FlatteRemain(false)
  ,_equalParameter(false)
{
  setUp(hypMap); 
}

Hyp3Lh::Hyp3Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps ) :
  Hyp2Lh(theLhPtr->getEventList(), hypMap, theStatesPtr, cacheAmps)
  ,_sigmaf980Hyp(false)
  ,_sigmaf1710Hyp(false)
  ,_sigmaf2200Hyp(false)
  ,_doHyp3(true)
  ,_nFitParams(0)
  ,_f980FlatteRemain(false)
  ,_equalParameter(false)
{
  setUp(hypMap); 
}

Hyp3Lh::~Hyp3Lh()
{;
}

bool  Hyp3Lh::equalChic0DecParams(){
  bool result=false; 
  bool equalRemainHyps=Hyp2Lh::equalChic0DecParams();
  if(!_doHyp3) return equalRemainHyps;
  _equalParameter=equalParams();
  DebugMsg << "equal parameter: "<< _equalParameter << endmsg;

  if(_equalParameter && equalRemainHyps) result=true;
  return result;
}

complex<double> Hyp3Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp2Lh::chi0DecAmps(theParamVal, theData);

  if(!_doHyp3) return result;

  if(_equalParameter){
    result+=_currentResultHyp3[_evtCounter];
    return result;
  }

  complex<double> currentResult(0.,0.);
 
  double sigmaMass=theParamVal.BwSigma.first;
  double sigmaWidth=theParamVal.BwSigma.second;

  if (_sigmaf980Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf980=theParamVal.ChiToSigmaf980;
 
    double f980_Mass=theParamVal.Flatf980;
    double f980_gPiPi=theParamVal.Flatf980gPiPi;
    double f980_gKK=theParamVal.Flatf980gKK;
    currentResult+= chiTof980_kf0_piAmp(theData, ChiToSigmaf980, f980_Mass, f980_gKK,  f980_gPiPi, sigmaMass, sigmaWidth);
  }

  if (_sigmaf1710Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710=theParamVal.ChiToSigmaf1710;
    double f1710Mass=theParamVal.Bwf1710.first;
    double f1710Width=theParamVal.Bwf1710.second;
    currentResult+=chiTof0_pif0_kAmp(theData, ChiToSigmaf1710, sigmaMass, sigmaWidth,  f1710Mass, f1710Width);
  }  

  if (_sigmaf2200Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf2200=theParamVal.ChiToSigmaf2200;
    double f2200_Mass=theParamVal.Bwf2200.first;
    double f2200_Width=theParamVal.Bwf2200.second;
    currentResult+=chiTof0_pif0_kAmp(theData, ChiToSigmaf2200, sigmaMass, sigmaWidth,  f2200_Mass, f2200_Width);
  }  

  if(_cacheAmps) _currentResultHyp3[_evtCounter]=currentResult; 
  result+=currentResult;
  return result;
}


void Hyp3Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp2Lh::setMnUsrParams(upar, startVal, errVal);
  if(!_doHyp3) return;


  if(_sigmaf980Hyp)  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::sigmaf980);
  if(_sigmaf1710Hyp) _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::sigmaf1710);
  if(_sigmaf2200Hyp) _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::sigmaf2200);

  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::sigma);

  if(_sigmaf2200Hyp && !_doHyp2)  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::f2200);   
}



int Hyp3Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }  

  int counter=Hyp2Lh::setFitParamVal(theParamVal, par);
  if(!_doHyp3) return counter;


  //sigma f1710    amplitude params
  if(_sigmaf980Hyp) counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::sigmaf980);
  if(_sigmaf1710Hyp) counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::sigmaf1710);
  if(_sigmaf2200Hyp) counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::sigmaf2200);
  
  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::sigma);  
  if(_sigmaf2200Hyp && !_doHyp2)  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::f2200);
     
  return counter;
}

unsigned int  Hyp3Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp2Lh::nFitParams();
  return _nFitParams+remainingFitParams;
}

void Hyp3Lh::print(std::ostream& os) const{
  os << "Hyp3Lh::print\n";
}

void Hyp3Lh::printCurrentFitResult(param2K2PiGam& theParamVal){
  Hyp2Lh::printCurrentFitResult(theParamVal);

  if(!_doHyp3) return;

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

void Hyp3Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }
  Hyp2Lh::dumpCurrentResult(os, theParamVal, suffix);

  if(!_doHyp3) return;

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

void Hyp3Lh::setUp(const std::map<const std::string, bool>& hypMap){

  std::string theKey="sigmaf980Hyp3";
  setHyps( hypMap, _sigmaf980Hyp, theKey);

  theKey="sigmaf1710Hyp3";
  setHyps( hypMap, _sigmaf1710Hyp, theKey);

  theKey="sigmaf2200Hyp3";
  setHyps( hypMap, _sigmaf2200Hyp, theKey);


  if (!_sigmaf980Hyp && !_sigmaf1710Hyp && !_sigmaf2200Hyp) _doHyp3=false;

  if (!_doHyp3) return;

  _massVec.push_back(paramEnum2K2PiGam::sigma); 

  if (_sigmaf980Hyp){ 
    _ampVec.push_back(paramEnum2K2PiGam::sigmaf980);
  }

  if (_sigmaf1710Hyp){ 
    _ampVec.push_back(paramEnum2K2PiGam::sigmaf1710);
  }

  if (_sigmaf2200Hyp){ 
    _ampVec.push_back(paramEnum2K2PiGam::sigmaf2200);
    if (!_doHyp2) _massVec.push_back(paramEnum2K2PiGam::f2200);  
  }

  // fill all other resonances

  if(_sigmaf980Hyp){
    _f980FlatteRemain=true;
  }
  
  if (_sigmaf1710Hyp){ 
    _massVecRemain.push_back(paramEnum2K2PiGam::f1710);
 }
  
  if (_sigmaf2200Hyp && _doHyp2){ 
    _massVecRemain.push_back(paramEnum2K2PiGam::f2200);
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



void Hyp3Lh::copyCurrentVals(Hyp3Lh* theLh){
  Hyp1Lh::copyCurrentVals(theLh);
  if(_cacheAmps){
    theLh->_currentResultHyp3=_currentResultHyp3;
  }  
}

bool Hyp3Lh::equalParams(){
  bool result=true;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  if (!compAmpParms( _ampVec )) return false;
  if (!compMassParms(_massVec)) return false;
  if (!compMassParms(_massVecRemain)) return false;
  if(_f980FlatteRemain){
    if (!compFlatteParms()) return false;
  }
  return result;
}
