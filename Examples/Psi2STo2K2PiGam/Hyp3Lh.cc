#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp3Lh::Hyp3Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap ) :
  Hyp2Lh(theEvtList, hypMap)
  ,_sigmaf980Hyp(true)
  ,_sigmaf1710Hyp(true)
  ,_sigmaf2200Hyp(true)  
  ,_disableHyp3(false)
  ,_nFitParams(0)
{
  setUp(hypMap); 
}

Hyp3Lh::Hyp3Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp2Lh(theLhPtr->getEventList(), hypMap)
  ,_sigmaf980Hyp(true)
  ,_sigmaf1710Hyp(true)
  ,_sigmaf2200Hyp(true)
  ,_disableHyp3(false)
  ,_nFitParams(0)
{
  setUp(hypMap); 
}

Hyp3Lh::~Hyp3Lh()
{;
}


complex<double> Hyp3Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp2Lh::chi0DecAmps(theParamVal, theData);

  if(_disableHyp3) return result;
 
  double sigmaMass=theParamVal.BwSigma.first;
  double sigmaWidth=theParamVal.BwSigma.second;

  if (_sigmaf980Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf980=theParamVal.ChiToSigmaf980;
 
    double f980_Mass=theParamVal.Flatf980;
    double f980_gPiPi=theParamVal.Flatf980gPiPi;
    double f980_gKK=theParamVal.Flatf980gKK;
    result+= chiTof980_kf0_piAmp(theData, ChiToSigmaf980, f980_Mass, f980_gKK,  f980_gPiPi, sigmaMass, sigmaWidth);
  }

  if (_sigmaf1710Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710=theParamVal.ChiToSigmaf1710;
    double f1710Mass=theParamVal.Bwf1710.first;
    double f1710Width=theParamVal.Bwf1710.second;
    result+=chiTof0_pif0_kAmp(theData, ChiToSigmaf1710, sigmaMass, sigmaWidth,  f1710Mass, f1710Width);
  }  

  if (_sigmaf2200Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf2200=theParamVal.ChiToSigmaf2200;
    double f2200_Mass=theParamVal.Bwf2200.first;
    double f2200_Width=theParamVal.Bwf2200.second;
    result+=chiTof0_pif0_kAmp(theData, ChiToSigmaf2200, sigmaMass, sigmaWidth,  f2200_Mass, f2200_Width);
  }  

  return result;
}


void Hyp3Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal){

  Hyp2Lh::setMnUsrParams(upar, startVal, errVal);
  if(_disableHyp3) return;


  if(_sigmaf980Hyp)  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::sigmaf980);
  if(_sigmaf1710Hyp) _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::sigmaf1710);
  if(_sigmaf2200Hyp) _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::sigmaf2200);

  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::sigma);

  if(_sigmaf2200Hyp && _disableHyp2)  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::f2200);   
}



int Hyp3Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }  

  int counter=Hyp2Lh::setFitParamVal(theParamVal, par);
  if(_disableHyp3) return counter;


  //sigma f1710    amplitude params
  if(_sigmaf980Hyp) counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::sigmaf980);
  if(_sigmaf1710Hyp) counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::sigmaf1710);
  if(_sigmaf2200Hyp) counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::sigmaf2200);
  
  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::sigma);  
  if(_sigmaf2200Hyp && _disableHyp2)  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::f2200);
     
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

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("sigmaf980Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf980Hyp= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf980Hyp <<endmsg;
  }
  else Alert << "hypothesis sigmaf980Hyp3 not set!!!" <<endmsg;

 iter= hypMap.find("sigmaf1710Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf1710Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf1710Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis sigmaf1710Hyp not set!!!" <<endmsg;

 iter= hypMap.find("sigmaf2200Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf2200Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf2200Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis sigmaf2200Hyp not set!!!" <<endmsg;

  if (!_sigmaf980Hyp && !_sigmaf1710Hyp && !_sigmaf2200Hyp) _disableHyp3=true;

  if (!_disableHyp3) _massVec.push_back(paramEnum2K2PiGam::sigma); 

  if (_sigmaf980Hyp){ 
    _ampVec.push_back(paramEnum2K2PiGam::sigmaf980);
  }

  if (_sigmaf1710Hyp){ 
    _ampVec.push_back(paramEnum2K2PiGam::sigmaf1710);
  }

  if (_sigmaf2200Hyp){ 
    _ampVec.push_back(paramEnum2K2PiGam::sigmaf2200);
    if (_disableHyp2) _massVec.push_back(paramEnum2K2PiGam::f2200);  
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
