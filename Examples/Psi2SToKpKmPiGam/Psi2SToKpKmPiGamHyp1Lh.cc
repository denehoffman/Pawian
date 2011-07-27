#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamHyp1Lh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamHyp1Lh::Psi2SToKpKmPiGamHyp1Lh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  Psi2SToKpKmPiGamBaseLh(theEvtList,  hypMap)
  ,_nFitParams(0)  
{
  _nFitParams=Psi2SToKpKmPiGamBaseLh::nFitParams();

  _ampVec.push_back(paramEnumKpKmPiGam::KappaK);

  _massVec.push_back(paramEnumKpKmPiGam::Kappa);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  std::vector<unsigned int>::iterator massIt; 
  for (massIt=_massVec.begin(); massIt!=_massVec.end(); ++massIt){
    _nFitParams+=2;
  }

}

Psi2SToKpKmPiGamHyp1Lh::Psi2SToKpKmPiGamHyp1Lh( boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Psi2SToKpKmPiGamBaseLh(theLhPtr->getEventList(), hypMap)
  ,_nFitParams(0)  
{
  _nFitParams=Psi2SToKpKmPiGamBaseLh::nFitParams();

  _ampVec.push_back(paramEnumKpKmPiGam::KappaK);

  _massVec.push_back(paramEnumKpKmPiGam::Kappa);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  std::vector<unsigned int>::iterator massIt; 
  for (massIt=_massVec.begin(); massIt!=_massVec.end(); ++massIt){
    _nFitParams+=2;
  }
}

Psi2SToKpKmPiGamHyp1Lh::~Psi2SToKpKmPiGamHyp1Lh()
{;
}


complex<double> Psi2SToKpKmPiGamHyp1Lh::calcDecAmp(complex<double>& inAmp,Spin lamChi, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){


  complex<double> result=Psi2SToKpKmPiGamBaseLh::calcDecAmp(inAmp, lamChi, theParamVal, theData);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKappaK=theParamVal.ChiToKappa_K;


  double KappaMass=theParamVal.BwKappa.first;
  double KappaWidth=theParamVal.BwKappa.second;

  //Chi to a0(980) pi0 and a0(980)->K+K-
  result+=K0_1400Amp(theData, ChiToKappaK, KappaMass, KappaWidth, lamChi);
  return result;
}



void Psi2SToKpKmPiGamHyp1Lh::setMnUsrParams(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal){

  Psi2SToKpKmPiGamBaseLh::setMnUsrParams(upar, startVal, errVal);

  _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal,paramEnumKpKmPiGam::KappaK);

  _fitParamsKpKmPiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::Kappa);
}


int Psi2SToKpKmPiGamHyp1Lh::setFitParamVal(paramKpKmPiGam& theParamVal, const std::vector<double>& par) {
  
  if (par.size()!= nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  int counter=0;

  counter= Psi2SToKpKmPiGamBaseLh::setFitParamVal(theParamVal, par);

  counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::KappaK);

  counter=_fitParamsKpKmPiGam.setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::Kappa);

  return counter;
}

unsigned int  Psi2SToKpKmPiGamHyp1Lh::nFitParams(){

  return _nFitParams;
}


void Psi2SToKpKmPiGamHyp1Lh::print(std::ostream& os) const{
  os << "Psi2SToKpKmPiGamHyp1Lh::print\n";
}

void Psi2SToKpKmPiGamHyp1Lh::printCurrentFitResult(paramKpKmPiGam& theParamVal){

  Psi2SToKpKmPiGamBaseLh::printCurrentFitResult(theParamVal);

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs = _fitParamsKpKmPiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< paramEnumKpKmPiGam::name(*itAmps) << endmsg;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParamsKpKmPiGam.ampMap(theParamVal, *itAmps);

      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }  
  }


  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    DebugMsg<< paramEnumKpKmPiGam::name(*itMasses) << endmsg;
    std::pair<double, double> tmpParam=_fitParamsKpKmPiGam.massPair(theParamVal, *itMasses);
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

}

void Psi2SToKpKmPiGamHyp1Lh::dumpCurrentResult(std::ostream& os, paramKpKmPiGam& theParamVal, std::string& suffix){

  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Psi2SToKpKmPiGamBaseLh::dumpCurrentResult(os, theParamVal, suffix);

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParamsKpKmPiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];

      std::string tmpStringDec=(*itJPCLS)->name()+paramEnumKpKmPiGam::name(*itAmps)+suffix;
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }  
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    std::string tmpStringMass=paramEnumKpKmPiGam::name(*itMasses)+suffix;

    std::pair<double, double> tmpParam=_fitParamsKpKmPiGam.massPair(theParamVal, *itMasses);
    os << tmpStringMass << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

}
