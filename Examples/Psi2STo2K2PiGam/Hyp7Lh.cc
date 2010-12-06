#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp7Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp7Lh::Hyp7Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp6Lh(theEvtList, theStates, hypMap )
{

}

Hyp7Lh::Hyp7Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp6Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
{

}

Hyp7Lh::~Hyp7Lh()
{;
}



complex<double> Hyp7Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp6Lh::chi0DecAmps(theParamVal, theData);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK892K1680=theParamVal.ChiToK892K1680;

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;
  double K_1_1680Mass=theParamVal.BwK_1_1680.first;
  double K_1_1680Width=theParamVal.BwK_1_1680.second;

  //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K


  result+=chiToK1K1Amp(theData, ChiToK892K1680, K892Mass, K892Width,  K_1_1680Mass, K_1_1680Width); 


  return result;
}




void Hyp7Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);
  Hyp6Lh::setMnUsrParams(upar, startVal, errVal);

    setMnUsrParamsDec(upar, startVal, errVal,"K892K_1_1680");
    setMnUsrParamsMass(upar, startVal, errVal, "K_1_1680");

}



int Hyp7Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
 
  int counter=Hyp6Lh::setFitParamVal(theParamVal, par);

  counter=setFitParamValDec(theParamVal, par, counter, "K892K_1_1680");

  counter=setFitParamValMass(theParamVal, par, counter, "K_1_1680");  

  return counter;
}

void Hyp7Lh::print(std::ostream& os) const{
  os << "Hyp7Lh::print\n";
}

void Hyp7Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp6Lh::printCurrentFitResult(theParamVal);

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK892K1680States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK892K1680States.begin(); itJPCLS!=JPCLSChiToK892K1680States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K892K_1_1680" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.ChiToK892K1680[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
      
    }
    
    DebugMsg<< "K_1_1680:" << endmsg;
    std::pair<double, double> tmpParamK_1_1680=theParamVal.BwK_1_1680;
    DebugMsg <<"\t mass:" << tmpParamK_1_1680.first <<"\t width:" << tmpParamK_1_1680.second  << endmsg;
}

void Hyp7Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp6Lh::dumpCurrentResult(os, theParamVal, suffix);

    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK892K1680States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK892K1680States.begin(); itJPCLS!=JPCLSChiToK892K1680States.end(); ++itJPCLS){
      std::string tmpStringDec=(*itJPCLS)->name()+"K892K_1_1680"+suffix;
      
      std::pair<double, double> tmpParam=theParamVal.ChiToK892K1680[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
    
    std::string tmpStringRes="K_1_1680"+suffix;
    std::pair<double, double> tmpParamK_1_1680=theParamVal.BwK_1_1680;
    os << tmpStringRes << "\t" << tmpParamK_1_1680.first <<"\t" << tmpParamK_1_1680.second  << std::endl;
}
