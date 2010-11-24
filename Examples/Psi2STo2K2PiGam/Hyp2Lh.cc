#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp2Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp2Lh::Hyp2Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, bool K1_1270Hyp, bool K0_1430_K0_1430Hyp) :
  Hyp1Lh(theEvtList,theStates, K1_1270Hyp, K0_1430_K0_1430Hyp)  
{
}

Hyp2Lh::Hyp2Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, bool K1_1270Hyp, bool K0_1430_K0_1430Hyp ) :
  Hyp1Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), K1_1270Hyp, K0_1430_K0_1430Hyp)
{
}

Hyp2Lh::~Hyp2Lh()
{;
}


complex<double> Hyp2Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f2200=theParamVal.ChiTof980f2200;
  double f2200Mass=theParamVal.Bwf2200.first;
  double f2200Width=theParamVal.Bwf2200.second;
  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  //Chi_c0 decay to f0(980) f0(2200) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980f2200Amp=chiTof980_pif0_kAmp(theData, ChiTof980f2200, f980_Mass, f980_gPiPi,  f980_gKK, f2200Mass, f2200Width);

  complex<double> result=Hyp1Lh::chi0DecAmps(theParamVal, theData)+ChiTof980f2200Amp;
                                                         
  return result;
}


void Hyp2Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  Hyp1Lh::setMnUsrParams(upar, startVal,  errVal);

  setMnUsrParamsDec(upar, startVal, errVal,"f980f2200");

  setMnUsrParamsMass(upar, startVal, errVal, "f2200");
}



int Hyp2Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
  if(_K0_1430_K0_1430Hyp) JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  
  if (par.size()< JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2
      +JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*8+15) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be more than" << JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2+JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*8+15 << endmsg;
    exit(1);
  }

  int counter=Hyp1Lh::setFitParamVal(theParamVal, par);

  //f980 f2200    amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "f980f2200");

  counter=setFitParamValMass(theParamVal, par, counter, "f2200");

  return counter;
}

void Hyp2Lh::print(std::ostream& os) const{
  os << "Hyp2Lh::print\n";
}

void Hyp2Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp1Lh::printCurrentFitResult(theParamVal);
//  print fit paramss
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "f980f2200" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980f2200[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

    DebugMsg<< "f2200_k:" << endmsg;
    std::pair<double, double> tmpParamf2200=theParamVal.Bwf2200;
    DebugMsg <<"\t mass:" << tmpParamf2200.first <<"\t width:" << tmpParamf2200.second  << endmsg;

}

void Hyp2Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  Hyp1Lh::dumpCurrentResult(os, theParamVal, suffix);

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    std::string tmpStringDec=(*itJPCLS)->name()+"f980f2200"+suffix;
    
    std::pair<double, double> tmpParam=theParamVal.ChiTof980f2200[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }
  
  std::string tmpStringRes="f2200mass"+suffix;
  std::pair<double, double> tmpParamf2200=theParamVal.Bwf2200;
  os << tmpStringRes << "\t" << tmpParamf2200.first <<"\t" << tmpParamf2200.second  << std::endl;
}
