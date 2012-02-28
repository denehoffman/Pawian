#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp2Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp2Lh::Hyp2Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, const std::map<const std::string, bool>& hypMap, bool chacheAmps ) :
  Hyp1Lh(theEvtList, hypMap, chacheAmps )
  ,_doHyp2(false)
  ,_nFitParams(0)
  ,_equalParameter(false)
{
  setUp(hypMap);
}

Hyp2Lh::Hyp2Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap, bool chacheAmps ) :
  Hyp1Lh(theLhPtr->getEventList(), hypMap, chacheAmps)
  ,_doHyp2(true)
  ,_nFitParams(0)
  ,_equalParameter(false)
{
  setUp(hypMap);
}

Hyp2Lh::~Hyp2Lh()
{
}

bool  Hyp2Lh::equalChic0DecParams(){
  bool result=false; 
  bool equalRemainHyps=Hyp1Lh::equalChic0DecParams();
  if(!_doHyp2) return equalRemainHyps;

  _equalParameter=equalParams();
  DebugMsg << "equal parameter: "<< _equalParameter << endmsg;

  if(_equalParameter && equalRemainHyps) result=true;
  return result;
}

complex<double> Hyp2Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  complex<double> result=Hyp1Lh::chi0DecAmps(theParamVal, theData);

  if(!_doHyp2) return result;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f2200=theParamVal.ChiTof980f2200;
  double f2200Mass=theParamVal.Bwf2200.first;
  double f2200Width=theParamVal.Bwf2200.second;
  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  //Chi_c0 decay to f0(980) f0(2200) -> (pi0 pi0) (K K) 
  result+=chiTof980_pif0_kAmp(theData, ChiTof980f2200, f980_Mass, f980_gPiPi,  f980_gKK, f2200Mass, f2200Width);

  return result;
}


void Hyp2Lh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal, param2K2PiGam& errVal){

  Hyp1Lh::setMnUsrParams(upar, startVal,  errVal);
  if(!_doHyp2) return;

  _fitParams2K2PiGam.setMnUsrParamsDec(upar, startVal, errVal, paramEnum2K2PiGam::f980f2200);
  _fitParams2K2PiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnum2K2PiGam::f2200);
}



int Hyp2Lh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }  

  int counter=Hyp1Lh::setFitParamVal(theParamVal, par);

  if(!_doHyp2) return counter;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  //f980 f2200    amplitude params
  counter=_fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::f980f2200);

  counter=_fitParams2K2PiGam.setFitParamValMass(theParamVal, par, counter, paramEnum2K2PiGam::f2200);

  return counter;
}


unsigned int  Hyp2Lh::nFitParams(){
  unsigned int remainingFitParams=Hyp1Lh::nFitParams(); 
  return _nFitParams+remainingFitParams;
}


void Hyp2Lh::print(std::ostream& os) const{
  os << "Hyp2Lh::print\n";
}

void Hyp2Lh::printCurrentFitResult(param2K2PiGam& theParamVal){

  Hyp1Lh::printCurrentFitResult(theParamVal);

  if(!_doHyp2) return;

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

void Hyp2Lh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp1Lh::dumpCurrentResult(os, theParamVal, suffix);

  if(!_doHyp2) return;

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


void Hyp2Lh::setUp(const std::map<const std::string, bool>& hypMap){

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("doHyp2");

  if (iter !=hypMap.end()){
    _doHyp2= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _doHyp2 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis doHyp2 not set!!!" <<endmsg; 

  if (!_doHyp2) return;

  _ampVec.push_back(paramEnum2K2PiGam::f980f2200);
  
  _massVec.push_back(paramEnum2K2PiGam::f2200);
  
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

bool Hyp2Lh::equalParams(){
  bool result=true;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  if (!compAmpParms( _ampVec )) return false;
  if (!compMassParms(_massVec)) return false;
  if (!compFlatteParms()) return false;

  return result;
}
