#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp4Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp4Lh::Hyp4Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp3Lh(theEvtList, theStates, hypMap )
{
}

Hyp4Lh::Hyp4Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap  ) :
  Hyp3Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
{
}

Hyp4Lh::~Hyp4Lh()
{;
}



complex<double> Hyp4Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1370_k=theParamVal.ChiTof980_pif1370_k;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1370_pi=theParamVal.ChiTof980_kf1370_pi;
  double f1370_kMass=theParamVal.Bwf1370_k.first;
  double f1370_kWidth=theParamVal.Bwf1370_k.second;
  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  //Chi_c0 decay to f0(980) f0(1370) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980_pif1370_kAmp=chiTof980_pif0_kAmp(theData, ChiTof980_pif1370_k, f980_Mass, f980_gPiPi,  f980_gKK, f1370_kMass, f1370_kWidth);
  complex<double> ChiTof980_kf1370_piAmp=chiTof980_kf0_piAmp(theData, ChiTof980_kf1370_pi, f980_Mass, f980_gKK,  f980_gPiPi, f1370_kMass, f1370_kWidth);

  complex<double> result=Hyp3Lh::chi0DecAmps(theParamVal, theData)+ChiTof980_pif1370_kAmp+ChiTof980_kf1370_piAmp;

  return result;
}


void Hyp4Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  Hyp3Lh::setMnUsrParams(upar, startVal, errVal);
  setMnUsrParamsDec(upar, startVal, errVal,"f980_pif1370_k");
  setMnUsrParamsDec(upar, startVal, errVal,"f980_kf1370_pi");
  setMnUsrParamsMass(upar, startVal, errVal, "f1370_k");  
}



int Hyp4Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();


  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
  if(_K0_1430_K0_1430Hyp) JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  
  if (par.size()< JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2
      +JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*14+19) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be more than" << JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2+JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*14+19 << endmsg;
    exit(1);
  }  

  int counter=Hyp3Lh::setFitParamVal(theParamVal, par);

  //sigma f1710    amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "f980_pif1370_k");
  counter=setFitParamValDec(theParamVal, par, counter, "f980_kf1370_pi");
  counter=setFitParamValMass(theParamVal, par, counter, "f1370_k");  

  return counter;
}

void Hyp4Lh::print(std::ostream& os) const{
  os << "Hyp4Lh::print\n";
}

void Hyp4Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp3Lh::printCurrentFitResult(theParamVal);

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "f980_pif1370_k" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_pif1370_k[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;

  }

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "f980_kf1370_pi" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_kf1370_pi[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;

  }

  DebugMsg<< "f1370_k:" << endmsg;
  std::pair<double, double> tmpParamf1370_k=theParamVal.Bwf1370_k;
  DebugMsg <<"\t mass:" << tmpParamf1370_k.first <<"\t width:" << tmpParamf1370_k.second  << endmsg;
}

void Hyp4Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }
  Hyp3Lh::dumpCurrentResult(os, theParamVal, suffix);

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    std::string tmpStringDec=(*itJPCLS)->name()+"f980_pif1370_k"+suffix;
    
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_pif1370_k[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    std::string tmpStringDec=(*itJPCLS)->name()+"f980_kf1370_pi"+suffix;
    
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_kf1370_pi[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

  
  std::string tmpStringRes="f1370_kmass"+suffix;
  std::pair<double, double> tmpParamf1370_k=theParamVal.Bwf1370_k;
  os << tmpStringRes << "\t" << tmpParamf1370_k.first <<"\t" << tmpParamf1370_k.second  << std::endl;
}
