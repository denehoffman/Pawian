#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp6Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp6Lh::Hyp6Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp5Lh(theEvtList, theStates, hypMap )
{

}

Hyp6Lh::Hyp6Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp5Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
{

}

Hyp6Lh::~Hyp6Lh()
{;
}



complex<double> Hyp6Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp5Lh::chi0DecAmps(theParamVal, theData);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1430K_0_1950=theParamVal.ChiToK_0_1430K_0_1950;


  double K_0_1430Mass=theParamVal.BwK_0_1430.first;
  double K_0_1430Width=theParamVal.BwK_0_1430.second;
  double K_0_1950Mass=theParamVal.BwK_0_1950.first;
  double K_0_1950Width=theParamVal.BwK_0_1950.second;

  //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K 
  complex<double> ChiToK1430K1950_Amp=chiTo2K_0_Amp(theData, ChiToK_0_1430K_0_1950, K_0_1430Mass, K_0_1430Width, K_0_1950Mass, K_0_1950Width); 


  result+=ChiToK1430K1950_Amp;

  return result;
}




void Hyp6Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);
  Hyp5Lh::setMnUsrParams(upar, startVal, errVal);

    setMnUsrParamsDec(upar, startVal, errVal,"K_0_1430K_0_1950");
    setMnUsrParamsMass(upar, startVal, errVal, "K_0_1950");

}



int Hyp6Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
 
  int counter=Hyp5Lh::setFitParamVal(theParamVal, par);

  counter=setFitParamValDec(theParamVal, par, counter, "K_0_1430K_0_1950");

  counter=setFitParamValMass(theParamVal, par, counter, "K_0_1950");  

  return counter;
}

void Hyp6Lh::print(std::ostream& os) const{
  os << "Hyp6Lh::print\n";
}

void Hyp6Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp5Lh::printCurrentFitResult(theParamVal);

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_0_1430K_0_1950" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.ChiToK_0_1430K_0_1950[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
      
    }
    
    DebugMsg<< "K_0_1950:" << endmsg;
    std::pair<double, double> tmpParamK_0_1950=theParamVal.BwK_0_1950;
    DebugMsg <<"\t mass:" << tmpParamK_0_1950.first <<"\t width:" << tmpParamK_0_1950.second  << endmsg;
}

void Hyp6Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp5Lh::dumpCurrentResult(os, theParamVal, suffix);

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      std::string tmpStringDec=(*itJPCLS)->name()+"K_0_1430K_0_1950"+suffix;
      
      std::pair<double, double> tmpParam=theParamVal.ChiToK_0_1430K_0_1950[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
    
    std::string tmpStringRes="K_0_1950"+suffix;
    std::pair<double, double> tmpParamK_0_1950=theParamVal.BwK_0_1950;
    os << tmpStringRes << "\t" << tmpParamK_0_1950.first <<"\t" << tmpParamK_0_1950.second  << std::endl;
}
