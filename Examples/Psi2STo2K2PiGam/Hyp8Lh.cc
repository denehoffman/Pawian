#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp8Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp8Lh::Hyp8Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp7Lh(theEvtList, theStates, hypMap )
{

}

Hyp8Lh::Hyp8Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp7Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
{

}

Hyp8Lh::~Hyp8Lh()
{;
}



complex<double> Hyp8Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp7Lh::chi0DecAmps(theParamVal, theData);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK892Pi=theParamVal.K_0_1460ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_0_1460ToK_0_1430Pi=theParamVal.K_0_1460ToK_0_1430Pi;

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;
  double K_0_1430Mass=theParamVal.BwK_0_1430.first;
  double K_0_1430Width=theParamVal.BwK_0_1430.second;
  double K_0_1460Mass=theParamVal.BwK_0_1460.first;
  double K_0_1460Width=theParamVal.BwK_0_1460.second;

  //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K


  result+=chiToK0KT0K1piKAmp(theData, K_0_1460ToK892Pi, K_0_1460Mass, K_0_1460Width, K892Mass, K892Width);
  result+=chiToK0KT0KpipiKAmp(theData, K_0_1460ToK_0_1430Pi, K_0_1460Mass, K_0_1460Width, K_0_1430Mass, K_0_1430Width); 


  return result;
}




void Hyp8Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);
  Hyp7Lh::setMnUsrParams(upar, startVal, errVal);

    setMnUsrParamsDec(upar, startVal, errVal,"K_0_1460ToK892Pi");
    setMnUsrParamsDec(upar, startVal, errVal, "K_0_1460ToK_0_1430Pi");
    setMnUsrParamsMass(upar, startVal, errVal, "K_0_1460");
}



int Hyp8Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
 
  int counter=Hyp7Lh::setFitParamVal(theParamVal, par);

  counter=setFitParamValDec(theParamVal, par, counter, "K_0_1460ToK892Pi");
  counter=setFitParamValDec(theParamVal, par, counter, "K_0_1460ToK_0_1430Pi");
  counter=setFitParamValMass(theParamVal, par, counter, "K_0_1460");  

  return counter;
}

void Hyp8Lh::print(std::ostream& os) const{
  os << "Hyp8Lh::print\n";
}

void Hyp8Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp7Lh::printCurrentFitResult(theParamVal);

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_0_1460ToK892Pi" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.K_0_1460ToK892Pi[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;

      DebugMsg<< (*itJPCLS)->name()<< "K_0_1460ToK_0_1430Pi" << endmsg;
      tmpParam=theParamVal.K_0_1460ToK_0_1430Pi[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
      
    }
    
    DebugMsg<< "K_0_1460:" << endmsg;
    std::pair<double, double> tmpParamK_0_1460=theParamVal.BwK_0_1460;
    DebugMsg <<"\t mass:" << tmpParamK_0_1460.first <<"\t width:" << tmpParamK_0_1460.second  << endmsg;
}

void Hyp8Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp7Lh::dumpCurrentResult(os, theParamVal, suffix);

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      std::string tmpStringDec=(*itJPCLS)->name()+"K_0_1460ToK892Pi"+suffix;
      
      std::pair<double, double> tmpParam=theParamVal.K_0_1460ToK892Pi[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }

    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      std::string tmpStringDec=(*itJPCLS)->name()+"K_0_1460ToK_0_1430Pi"+suffix;
      
      std::pair<double, double> tmpParam=theParamVal.K_0_1460ToK_0_1430Pi[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }

    
    std::string tmpStringRes="K_0_1460"+suffix;
    std::pair<double, double> tmpParamK_0_1460=theParamVal.BwK_0_1460;
    os << tmpStringRes << "\t" << tmpParamK_0_1460.first <<"\t" << tmpParamK_0_1460.second  << std::endl;
}
