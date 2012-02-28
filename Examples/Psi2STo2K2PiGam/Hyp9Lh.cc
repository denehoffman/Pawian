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
  ,_Pi1800Tof980PiHyp9(true)
  ,_Pi1800Tof1370PiHyp9(true)
  ,_Pi1800ToKappaKHyp9(true)
  ,_Pi1800Pi0ToK892KHyp9(true)
  ,_Pi3000Pi0ToK892KHyp9(true)
  ,_Pi3000Pi0ToK_0_1950KHyp9(true)
  ,_Pi_2_2285Tof1710PiHyp9(true)
  ,_Pi_2_2285ToK892KHyp9(true)
  ,_Pi_2_2285ToK_0_1430KHyp9(true)
  ,_Pi_2_2285ToK_2_1430KHyp9(true)
  ,_f980f_2_2300Hyp9(true)
  ,_f_2_2300sigmaHyp9(true)
  ,_K_2_1770ToK_2_1430PiHyp9(true)
  ,_doHyp9(true)
  ,_nFitParams(0)
  ,_f980FlatteRemain(false)
  ,_evtCounter(0)
  ,_equalParameter(false)
  ,_equalPi_2_2285ToK_2_1430Params(false)
  ,_equalPi_2_2285ToK892KParams(false)
{
  setUp(hypMap); 
}

Hyp9Lh::Hyp9Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp8Lh(theLhPtr->getEventList(), hypMap)
  ,_Pi_2_1670Tof_2_1270PiHyp9(true)
  ,_Pi_2_1670ToK892KHyp9(true)
  ,_Pi1800Tof980PiHyp9(true)
  ,_Pi1800Tof1370PiHyp9(true)
  ,_Pi1800ToKappaKHyp9(true)
  ,_Pi1800Pi0ToK892KHyp9(true)
  ,_Pi3000Pi0ToK892KHyp9(true)
  ,_Pi3000Pi0ToK_0_1950KHyp9(true)
  ,_Pi_2_2285Tof1710PiHyp9(true)
  ,_Pi_2_2285ToK892KHyp9(true)
  ,_Pi_2_2285ToK_0_1430KHyp9(true)
  ,_Pi_2_2285ToK_2_1430KHyp9(true)
  ,_f980f_2_2300Hyp9(true)
  ,_f_2_2300sigmaHyp9(true)
  ,_K_2_1770ToK_2_1430PiHyp9(true)
  ,_doHyp9(true)
  ,_nFitParams(0)
  ,_f980FlatteRemain(false)
  ,_evtCounter(0)
  ,_equalParameter(false)
  ,_equalPi_2_2285ToK_2_1430Params(false)
  ,_equalPi_2_2285ToK892KParams(false)
{
  setUp(hypMap); 
}

Hyp9Lh::~Hyp9Lh()
{;
}

bool  Hyp9Lh::equalChic0DecParams(){
  bool result=false; 
  bool equalRemainHyps=Hyp8Lh::equalChic0DecParams();
  if (!_doHyp9) return equalRemainHyps;
  _equalParameter=equalParams();
  DebugMsg << "equal parameter: "<< _equalParameter << endmsg;
  if (compAmpParms(_ampVecPi_2_2285ToK_2_1430) && compMassParms(_massVecPi_2_2285ToK_2_1430)) _equalPi_2_2285ToK_2_1430Params=true;
  if (compAmpParms(_ampVecPi_2_2285ToK892K) && compMassParms(_massVecPi_2_2285ToK892K)) _equalPi_2_2285ToK892KParams=true;
  if(_equalParameter && equalRemainHyps) result=true;
  return result;
}

complex<double> Hyp9Lh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp8Lh::chi0DecAmps(theParamVal, theData);

  if (!_doHyp9) return result;

  if (_evtCounter==0){
    _equalParameter=equalParams();
    if (compAmpParms(_ampVecPi_2_2285ToK_2_1430) && compMassParms(_massVecPi_2_2285ToK_2_1430)) _equalPi_2_2285ToK_2_1430Params=true;
    if (compAmpParms(_ampVecPi_2_2285ToK892K) && compMassParms(_massVecPi_2_2285ToK892K)) _equalPi_2_2285ToK892KParams=true;
    DebugMsg << "equal parameter: "<< _equalParameter << endmsg;

  } 

  if(_equalParameter){
    result+=_currentResultHyp9[_evtCounter];
    _evtCounter++;
    return result;
  }

  complex<double> currentResult(0.,0.);

  if(_Pi_2_1670Tof_2_1270PiHyp9 || _Pi_2_1670ToK892KHyp9){
    double Pi_2_1670Mass=theParamVal.BwPi_2_1670.first;
    double Pi_2_1670Width=theParamVal.BwPi_2_1670.second;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi_2_1670Pi=theParamVal.ChiToPi_2_1670Pi; 
    
    if(_Pi_2_1670Tof_2_1270PiHyp9){
      double f_2_1270Mass=theParamVal.Bwf_2_1270.first;
      double f_2_1270Width=theParamVal.Bwf_2_1270.second;
      
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_1670Tof_2_1270Pi=theParamVal.Pi_2_1670Tof_2_1270Pi;
      
      currentResult+=chiToPi2Pi0Tof2PiAmp(theData, ChiToPi_2_1670Pi, Pi_2_1670Tof_2_1270Pi, Pi_2_1670Mass, Pi_2_1670Width, f_2_1270Mass, f_2_1270Width);
      
    }
    
    if(_Pi_2_1670ToK892KHyp9){
      double K892Mass=theParamVal.BwK892.first;
      double K892Width=theParamVal.BwK892.second;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_1670ToK892K=theParamVal.Pi_2_1670ToK892K;
      
      currentResult+=chiToPi2Pi0ToKstarKAmp(theData, ChiToPi_2_1670Pi, Pi_2_1670ToK892K, Pi_2_1670Mass, Pi_2_1670Width, K892Mass, K892Width);
      
    }
  }

  if(_Pi1800Tof980PiHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0Tof980=theParamVal.ChiToPi1800Pi0Tof980;
    double Pi1800Mass=theParamVal.BwPi1800.first;
    double Pi1800Width=theParamVal.BwPi1800.second;
    double f980Mass=theParamVal.Flatf980;
    double f980gPiPi=theParamVal.Flatf980gPiPi;
    double f980gKK=theParamVal.Flatf980gKK;

    currentResult+=chiToPi0Pi0Tof980kkAmp(theData, ChiToPi1800Pi0Tof980, f980Mass, f980gKK, f980gPiPi, Pi1800Mass, Pi1800Width);
  }

  if(_Pi1800Tof1370PiHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0Tof1370=theParamVal.ChiToPi1800Pi0Tof1370;
    double Pi1800Mass=theParamVal.BwPi1800.first;
    double Pi1800Width=theParamVal.BwPi1800.second;
    double f1360Mass=theParamVal.Bwf1360.first;
    double f1360Width=theParamVal.Bwf1360.second;
    
    currentResult+=chiToPi0Pi0Tof0Amp(theData, ChiToPi1800Pi0Tof1370, Pi1800Mass, Pi1800Width, f1360Mass, f1360Width);
  }

  if(_Pi1800ToKappaKHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0ToKappa=theParamVal.ChiToPi1800Pi0ToKappa;
    double Pi1800Mass=theParamVal.BwPi1800.first;
    double Pi1800Width=theParamVal.BwPi1800.second;
    double KappaMass=theParamVal.BwKappa.first;
    double KappaWidth=theParamVal.BwKappa.second;
    
    currentResult+=chiToPi0Pi0ToK0KAmp(theData, ChiToPi1800Pi0ToKappa, Pi1800Mass, Pi1800Width, KappaMass, KappaWidth);
  }

  if(_Pi1800Pi0ToK892KHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi1800Pi0ToK892K=theParamVal.ChiToPi1800Pi0ToK892K;
    double Pi1800Mass=theParamVal.BwPi1800.first;
    double Pi1800Width=theParamVal.BwPi1800.second;
    double K892Mass=theParamVal.BwK892.first;
    double K892Width=theParamVal.BwK892.second;
    
    currentResult+=chiToPi0PiToKstarKAmp(theData, ChiToPi1800Pi0ToK892K, Pi1800Mass, Pi1800Width, K892Mass, K892Width);

  }

  if(_Pi3000Pi0ToK892KHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi3000Pi0ToK892K=theParamVal.ChiToPi3000Pi0ToK892K;
    double Pi3000Mass=theParamVal.BwPi3000.first;
    double Pi3000Width=theParamVal.BwPi3000.second;
    double K892Mass=theParamVal.BwK892.first;
    double K892Width=theParamVal.BwK892.second;
    
    currentResult+=chiToPi0PiToKstarKAmp(theData, ChiToPi3000Pi0ToK892K, Pi3000Mass, Pi3000Width, K892Mass, K892Width);

  }

  if(_Pi3000Pi0ToK_0_1950KHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi3000Pi0ToK_0_1950K=theParamVal.ChiToPi3000Pi0ToK_0_1950K;
    double Pi3000Mass=theParamVal.BwPi3000.first;
    double Pi3000Width=theParamVal.BwPi3000.second;
    double K_0_1950Mass=theParamVal.BwK_0_1950.first;
    double K_0_1950Width=theParamVal.BwK_0_1950.second;

    currentResult+=chiToPi0Pi0ToK0KAmp(theData, ChiToPi3000Pi0ToK_0_1950K, Pi3000Mass, Pi3000Width, K_0_1950Mass, K_0_1950Width);
  }

  if (_Pi_2_2285Tof1710PiHyp9 || _Pi_2_2285ToK892KHyp9 || _Pi_2_2285ToK_0_1430KHyp9 || _Pi_2_2285ToK_2_1430KHyp9){
    double Pi_2_2285Mass=theParamVal.BwPi_2_2285.first;
    double Pi_2_2285Width=theParamVal.BwPi_2_2285.second;
std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToPi_2_2285Pi=theParamVal.ChiToPi_2_2285Pi;
    if(_Pi_2_2285Tof1710PiHyp9){

    double f1710Mass=theParamVal.Bwf1710.first;
    double f1710Width=theParamVal.Bwf1710.second;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285Tof1700Pi=theParamVal.Pi_2_2285Tof1700Pi;

    currentResult+=chiToPi2Pi0Tof0PiAmp(theData, ChiToPi_2_2285Pi, Pi_2_2285Tof1700Pi, Pi_2_2285Mass, Pi_2_2285Width, f1710Mass, f1710Width);
    }

    if(_Pi_2_2285ToK892KHyp9){
      if(_equalPi_2_2285ToK892KParams) currentResult+=_currentResultPi_2_2285ToK892KHyp9[_evtCounter];
      else{
	std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285ToK892KPi=theParamVal.Pi_2_2285ToK892KPi;
	double K892Mass=theParamVal.BwK892.first;
	double K892Width=theParamVal.BwK892.second;
	complex<double> tmpPi_2_2285ToK892KAmp=chiToPi2Pi0ToKstarKAmp(theData, ChiToPi_2_2285Pi, Pi_2_2285ToK892KPi, Pi_2_2285Mass, Pi_2_2285Width, K892Mass, K892Width, 1);
	_currentResultPi_2_2285ToK892KHyp9[_evtCounter]=tmpPi_2_2285ToK892KAmp;	
	currentResult+=tmpPi_2_2285ToK892KAmp;
      }
    }

    if(_Pi_2_2285ToK_0_1430KHyp9){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285ToK_0_1430K=theParamVal.Pi_2_2285ToK_0_1430K;
    double K_0_1430Mass=theParamVal.BwK_0_1430.first;
    double K_0_1430Width=theParamVal.BwK_0_1430.second;

    currentResult+=chiToPi2Pi0ToKstarKAmp(theData, ChiToPi_2_2285Pi, Pi_2_2285ToK_0_1430K, Pi_2_2285Mass, Pi_2_2285Width, K_0_1430Mass, K_0_1430Width,0 );
    }

    if(_Pi_2_2285ToK_2_1430KHyp9){
      
      if(_equalPi_2_2285ToK_2_1430Params)  currentResult+=_currentResultPi_2_2285ToK_2_1430Hyp9[_evtCounter];
      else{
	std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > Pi_2_2285ToK_2_1430K=theParamVal.Pi_2_2285ToK_2_1430K;
	double K_2_1430Mass=theParamVal.BwK_2_1430.first;
	double K_2_1430Width=theParamVal.BwK_2_1430.second;
	
	complex<double> tmpPi_2_2285ToK_2_1430Amp=chiToPi2Pi0ToKstarKAmp(theData, ChiToPi_2_2285Pi, Pi_2_2285ToK_2_1430K, Pi_2_2285Mass, Pi_2_2285Width, K_2_1430Mass, K_2_1430Width,2 );
	_currentResultPi_2_2285ToK_2_1430Hyp9[_evtCounter]=tmpPi_2_2285ToK_2_1430Amp;
	currentResult+=tmpPi_2_2285ToK_2_1430Amp;
      }
    }
  }

  if (_f980f_2_2300Hyp9 || _f_2_2300sigmaHyp9){
    double f_2_2300Mass=theParamVal.Bwf_2_2300.first;
    double f_2_2300Width=theParamVal.Bwf_2_2300.second;

    if (_f980f_2_2300Hyp9){
    double f980Mass=theParamVal.Flatf980;
    double f980gPiPi=theParamVal.Flatf980gPiPi;
    double f980gKK=theParamVal.Flatf980gKK;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f_2_2300=theParamVal.ChiTof980f_2_2300;

    currentResult+=chiTof980_pif2_kAmp(theData, ChiTof980f_2_2300, f980Mass, f980gPiPi, f980gKK, f_2_2300Mass, f_2_2300Width);
    }

    if (_f_2_2300sigmaHyp9){
    double sigmaMass=theParamVal.BwSigma.first;
    double sigmaWidth=theParamVal.BwSigma.second;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof_2_2300sigma=theParamVal.ChiTof_2_2300sigma;

    currentResult+=chiTof0_pif2_kAmp(theData, ChiTof_2_2300sigma, sigmaMass, sigmaWidth, f_2_2300Mass, f_2_2300Width);
    }
  }



  if (_K_2_1770ToK_2_1430PiHyp9){
    double K_2_1770Mass=theParamVal.BwK_2_1770.first;
    double K_2_1770Width=theParamVal.BwK_2_1770.second;
    double K_2_1430Mass=theParamVal.BwK_2_1430.first;
    double K_2_1430Width=theParamVal.BwK_2_1430.second;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_2_1770K=theParamVal.ChiToK_2_1770K;
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_2_1770ToK_2_1430Pi=theParamVal.K_2_1770ToK_2_1430Pi;

    currentResult+=chiToK2KToK2PiAmp(theData, ChiToK_2_1770K, K_2_1770ToK_2_1430Pi, K_2_1770Mass, K_2_1770Width, K_2_1430Mass, K_2_1430Width);
  }

  _currentResultHyp9[_evtCounter]=currentResult; 
  _evtCounter++;
  result+=currentResult;
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
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi_2_1670Tof_2_1270PiHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_1670Tof_2_1270PiHyp9 not set!!!" <<endmsg;
    exit(0);
  }

  iter= hypMap.find("Pi_2_1670ToK892KHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_1670ToK892KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi_2_1670ToK892KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_1670ToK892KHyp9 not set!!!" <<endmsg;
    exit(0);
  }  

  iter= hypMap.find("Pi1800Tof980PiHyp9");
  
  if (iter !=hypMap.end()){
    _Pi1800Tof980PiHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi1800Tof980PiHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi1800Tof980PiHyp9 not set!!!" <<endmsg;
    exit(0);
  }  

  iter= hypMap.find("Pi1800ToKappaKHyp9");
    if (iter !=hypMap.end()){
    _Pi1800ToKappaKHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi1800ToKappaKHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi1800ToKappaKHyp9 not set!!!" <<endmsg;
    exit(0);
  } 


  iter= hypMap.find("Pi1800Tof1370PiHyp9");
  if (iter !=hypMap.end()){
    _Pi1800Tof1370PiHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi1800Tof1370PiHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi1800Tof1370PiHyp9 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("Pi1800Pi0ToK892KHyp9");
  if (iter !=hypMap.end()){
    _Pi1800Pi0ToK892KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi1800Pi0ToK892KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi1800Pi0ToK892KHyp9 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("Pi3000Pi0ToK892KHyp9");
  if (iter !=hypMap.end()){
    _Pi3000Pi0ToK892KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi3000Pi0ToK892KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi30000Pi0ToK892KHyp9 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("Pi3000Pi0ToK_0_1950KHyp9");
  if (iter !=hypMap.end()){
    _Pi3000Pi0ToK_0_1950KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi3000Pi0ToK892KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi3000Pi0ToK_0_1950KHyp9 not set!!!" <<endmsg;
    exit(0);
  } 

  iter= hypMap.find("Pi_2_2285Tof1710PiHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_2285Tof1710PiHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi_2_2285Tof1710PiHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_2285Tof1710PiHyp9 not set!!!" <<endmsg;
    exit(0);
  }


  iter= hypMap.find("Pi_2_2285ToK892KHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_2285ToK892KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi_2_2285ToK892KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_2285ToK892KHyp9 not set!!!" <<endmsg;
    exit(0);
  }

  iter= hypMap.find("Pi_2_2285ToK_0_1430KHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_2285ToK_0_1430KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi_2_2285ToK_0_1430KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_2285ToK_0_1430KHyp9 not set!!!" <<endmsg;
    exit(0);
  }




  iter= hypMap.find("Pi_2_2285ToK_2_1430KHyp9");
  
  if (iter !=hypMap.end()){
    _Pi_2_2285ToK_2_1430KHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _Pi_2_2285ToK_2_1430KHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "Pi_2_2285ToK_2_1430KHyp9 not set!!!" <<endmsg;
    exit(0);
  }



  iter= hypMap.find("f980f_2_2300Hyp9");
  
  if (iter !=hypMap.end()){
    _f980f_2_2300Hyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _f980f_2_2300Hyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "f980f_2_2300Hyp9 not set!!!" <<endmsg;
    exit(0);
  }

  iter= hypMap.find("f_2_2300sigmaHyp9");
  
  if (iter !=hypMap.end()){
    _f_2_2300sigmaHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _f_2_2300sigmaHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "f_2_2300sigmaHyp9 not set!!!" <<endmsg;
    exit(0);
  }

  iter= hypMap.find("K_2_1770ToK_2_1430PiHyp9");
  
  if (iter !=hypMap.end()){
    _K_2_1770ToK_2_1430PiHyp9= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << _K_2_1770ToK_2_1430PiHyp9 <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << "K_2_1770ToK_2_1430PiHyp9 not set!!!" <<endmsg;
    exit(0);
  }


  if (_Pi_2_1670Tof_2_1270PiHyp9 || _Pi_2_1670ToK892KHyp9){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToPi_2_1670Pi);
    _massVec.push_back(paramEnum2K2PiGam::Pi_2_1670);
    
    if (_Pi_2_1670Tof_2_1270PiHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::Pi_2_1670Tof_2_1270Pi);
      if (!_f980f_2_1270Hyp4){
	_massVec.push_back(paramEnum2K2PiGam::f_2_1270);
      }
    }
    
    if (_Pi_2_1670ToK892KHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::Pi_2_1670ToK892K);
    }
  }


  if (_Pi1800Tof980PiHyp9 || _Pi1800Tof1370PiHyp9 || _Pi1800ToKappaKHyp9 || _Pi1800ToKappaKHyp9){

    _massVec.push_back(paramEnum2K2PiGam::Pi1800);

    if (_Pi1800Tof980PiHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToPi1800Pi0Tof980);
    }
    
    if (_Pi1800Tof1370PiHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToPi1800Pi0Tof1370);
      _massVec.push_back(paramEnum2K2PiGam::f1360);
    }
    
    
    if(_Pi1800ToKappaKHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToPi1800Pi0ToKappa);
    }

    if (_Pi1800Pi0ToK892KHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::ChiToPi1800Pi0ToK892K);
    }

  }

    if (_Pi3000Pi0ToK892KHyp9 || _Pi3000Pi0ToK_0_1950KHyp9){
      _massVec.push_back(paramEnum2K2PiGam::Pi3000);

      if (_Pi3000Pi0ToK892KHyp9){
	_ampVec.push_back(paramEnum2K2PiGam::ChiToPi3000Pi0ToK892K);
      }
      if (_Pi3000Pi0ToK_0_1950KHyp9){
	_ampVec.push_back(paramEnum2K2PiGam::ChiToPi3000Pi0ToK_0_1950K);
	if(!_K_0_1430K_0_1950Hyp6 && !_KappaK_0_1950Hyp6){
	  _massVec.push_back(paramEnum2K2PiGam::K_0_1950);
	}
      }
    }

  if (_Pi_2_2285Tof1710PiHyp9 || _Pi_2_2285ToK892KHyp9 || _Pi_2_2285ToK_0_1430KHyp9 || _Pi_2_2285ToK_2_1430KHyp9){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToPi_2_2285Pi);
    _massVec.push_back(paramEnum2K2PiGam::Pi_2_2285);

    if (_Pi_2_2285Tof1710PiHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::Pi_2_2285Tof1700Pi);
    }

    if (_Pi_2_2285ToK892KHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::Pi_2_2285ToK892KPi);
      _ampVecPi_2_2285ToK892K.push_back(paramEnum2K2PiGam::ChiToPi_2_2285Pi);
      _ampVecPi_2_2285ToK892K.push_back(paramEnum2K2PiGam::Pi_2_2285ToK892KPi);
      _massVecPi_2_2285ToK892K.push_back(paramEnum2K2PiGam::Pi_2_2285);
      _massVecPi_2_2285ToK892K.push_back(paramEnum2K2PiGam::K892);
    }

    if(_Pi_2_2285ToK_0_1430KHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::Pi_2_2285ToK_0_1430K);
    }

    if(_Pi_2_2285ToK_2_1430KHyp9){
      _ampVec.push_back(paramEnum2K2PiGam::Pi_2_2285ToK_2_1430K);
      if(!_K2_1430_K2_1430Hyp && !_K0_1430_K2_1430Hyp) _massVec.push_back(paramEnum2K2PiGam::K_2_1430);

      _ampVecPi_2_2285ToK_2_1430.push_back(paramEnum2K2PiGam::ChiToPi_2_2285Pi);
      _ampVecPi_2_2285ToK_2_1430.push_back(paramEnum2K2PiGam::Pi_2_2285ToK_2_1430K);
      _massVecPi_2_2285ToK_2_1430.push_back(paramEnum2K2PiGam::Pi_2_2285);
      _massVecPi_2_2285ToK_2_1430.push_back(paramEnum2K2PiGam::K_2_1430);
    }
  }

  if (_f980f_2_2300Hyp9 || _f_2_2300sigmaHyp9){

    _massVec.push_back(paramEnum2K2PiGam::f_2_2300);
    if(_f980f_2_2300Hyp9){
    _ampVec.push_back(paramEnum2K2PiGam::ChiTof980f_2_2300);
    }
    if(_f_2_2300sigmaHyp9){
    _ampVec.push_back(paramEnum2K2PiGam::ChiTof_2_2300sigma);
    }
  }

  if (_K_2_1770ToK_2_1430PiHyp9){
    _ampVec.push_back(paramEnum2K2PiGam::ChiToK_2_1770K);
    _ampVec.push_back(paramEnum2K2PiGam::K_2_1770ToK_2_1430Pi);
    _massVec.push_back(paramEnum2K2PiGam::K_2_1770);
  }


  // fill all other resonances
  if (_Pi_2_2285Tof1710PiHyp9){
      _massVecRemain.push_back(paramEnum2K2PiGam::f1710);
  }

  if (_Pi_2_1670Tof_2_1270PiHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::f_2_1270);
  }

  if (_Pi_2_1670ToK892KHyp9 || _Pi1800Pi0ToK892KHyp9 || _Pi3000Pi0ToK892KHyp9 || _Pi_2_2285ToK892KHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::K892);
  }

  if (_Pi_2_2285ToK_0_1430KHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::K_0_1430);
  }

  if (_Pi_2_2285ToK_2_1430KHyp9 || _K_2_1770ToK_2_1430PiHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::K_2_1430);
  }

  if (_f_2_2300sigmaHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::sigma);
  }

  if(_Pi1800ToKappaKHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::Kappa);
  }

  if(_Pi3000Pi0ToK_0_1950KHyp9){
    _massVecRemain.push_back(paramEnum2K2PiGam::K_0_1950);
  }

  if(_Pi1800Tof980PiHyp9 || _f980f_2_2300Hyp9){
    _f980FlatteRemain=true;
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


void Hyp9Lh::copyCurrentVals(Hyp9Lh* theLh){
  Hyp8Lh::copyCurrentVals(theLh);
  std::map<unsigned int, complex<double> > newResult; 
  std::map<unsigned int, complex<double> >::iterator it;
  for (it= _currentResultHyp9.begin(); it!= _currentResultHyp9.end(); ++it){
    newResult[it->first]=it->second;
  } 

 std::map<unsigned int, complex<double> > newResultPi_2_2285ToK_2_1430; 
  for (it= _currentResultPi_2_2285ToK_2_1430Hyp9.begin(); it!= _currentResultPi_2_2285ToK_2_1430Hyp9.end(); ++it){
    newResultPi_2_2285ToK_2_1430[it->first]=it->second;
  }

 std::map<unsigned int, complex<double> > newResultPi_2_2285ToK892K; 
  for (it= _currentResultPi_2_2285ToK892KHyp9.begin(); it!= _currentResultPi_2_2285ToK892KHyp9.end(); ++it){
    newResultPi_2_2285ToK892K[it->first]=it->second;
  }

  theLh->_currentResultHyp9=newResult;
  theLh->_currentResultPi_2_2285ToK_2_1430Hyp9=newResultPi_2_2285ToK_2_1430; 
  theLh->_currentResultPi_2_2285ToK892KHyp9=newResultPi_2_2285ToK892K; 
}

bool Hyp9Lh::equalParams(){
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
