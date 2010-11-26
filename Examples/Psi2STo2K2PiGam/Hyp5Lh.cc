#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp5Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"


Hyp5Lh::Hyp5Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp4Lh(theEvtList, theStates, hypMap )
  , _disableHyp5(false)
{

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("disableHyp5");

  if (iter !=hypMap.end()){
    _disableHyp5= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _disableHyp5 <<endmsg;
  }
  else Alert << "hypothesis disableHyp5 not set!!!" <<endmsg;
}

Hyp5Lh::Hyp5Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp4Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
  , _disableHyp5(false)
{

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("disableHyp5");

  if (iter !=hypMap.end()){
    _disableHyp5= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _disableHyp5 <<endmsg;
  }
  else Alert << "hypothesis " << iter->first << "not set!!!" <<endmsg;
}

Hyp5Lh::~Hyp5Lh()
{;
}



complex<double> Hyp5Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp4Lh::chi0DecAmps(theParamVal, theData);
  if (_disableHyp5) return result;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_2400ToKf980=theParamVal.ChiToK_0_2400ToKf980;
  double K_0_2400Mass=theParamVal.BwK_0_2400.first;
  double K_0_2400Width=theParamVal.BwK_0_2400.second;
  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  //Chi_c0 decay to K_0_2400 K  -> (K f0(980) ) K -> (K pi0 pi0 ) K 
  complex<double> ChiToK2400KToKf980K_Amp=chiToK_0_KToKf980KAmp(theData, ChiToK_0_2400ToKf980, f980_Mass, f980_gPiPi,  f980_gKK, K_0_2400Mass, K_0_2400Width); 


  result+=ChiToK2400KToKf980K_Amp;

  return result;
}




void Hyp5Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);
    Hyp4Lh::setMnUsrParams(upar, startVal, errVal);

  if (!_disableHyp5){
    setMnUsrParamsDec(upar, startVal, errVal,"K_0_2400KToKf980K");
    setMnUsrParamsMass(upar, startVal, errVal, "K_0_2400");
  }  
}



int Hyp5Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
 
  int counter=Hyp4Lh::setFitParamVal(theParamVal, par);

  if (_disableHyp5) return counter;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();


  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
  if(_K0_1430_K0_1430Hyp) JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
  
  if (par.size()< JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2
      +JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiToK0K0States.size()*2+JPCLSChiTof0f0States.size()*12+21) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be more than" << JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2+JPCLSK1400ToKst1PiStates.size()*2
      +JPCLSChiToK0K0States.size()*2+JPCLSChiTof0f0States.size()*12+21 << endmsg;
    exit(1);
  } 


  counter=setFitParamValDec(theParamVal, par, counter, "K_0_2400KToKf980K");

  counter=setFitParamValMass(theParamVal, par, counter, "K_0_2400");  

  return counter;
}

void Hyp5Lh::print(std::ostream& os) const{
  os << "Hyp5Lh::print\n";
}

void Hyp5Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp4Lh::printCurrentFitResult(theParamVal);

  if (!_disableHyp5){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_0_2400KToKf980K" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.ChiToK_0_2400ToKf980[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
      
    }
    
    DebugMsg<< "K_0_2400:" << endmsg;
    std::pair<double, double> tmpParamK_0_2400=theParamVal.BwK_0_2400;
    DebugMsg <<"\t mass:" << tmpParamK_0_2400.first <<"\t width:" << tmpParamK_0_2400.second  << endmsg;
  }
}

void Hyp5Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp4Lh::dumpCurrentResult(os, theParamVal, suffix);

  if (!_disableHyp5){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiToK0K0States.begin(); itJPCLS!=JPCLSChiToK0K0States.end(); ++itJPCLS){
      std::string tmpStringDec=(*itJPCLS)->name()+"K_0_2400KToKf980K"+suffix;
      
      std::pair<double, double> tmpParam=theParamVal.ChiToK_0_2400ToKf980[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
    
    std::string tmpStringRes="K_0_2400"+suffix;
    std::pair<double, double> tmpParamK_0_2400=theParamVal.BwK_0_2400;
    os << tmpStringRes << "\t" << tmpParamK_0_2400.first <<"\t" << tmpParamK_0_2400.second  << std::endl;
  }
}
