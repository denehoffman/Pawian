#include <getopt.h>
#include <fstream>

#include "Examples/Psi2STo2K2PiGam/Hyp1Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp1Lh::Hyp1Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  AbsPsi2STo2K2PiGamLh(theEvtList,theStates)
  ,_K1_1270Hyp(true)
  ,_K1_1400Hyp(false)
  ,_K0_1430_K0_1430Hyp(true)
{
 
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K1_1270Hyp");

  if (iter !=hypMap.end()){
    _K1_1270Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1270Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1270Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K1_1400Hyp");
  if (iter !=hypMap.end()){
    _K1_1400Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1400Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1400Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K0_1430_K0_1430Hyp");
  if (iter !=hypMap.end()){
    _K0_1430_K0_1430Hyp= iter->second;
    _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1270Hyp <<endmsg;
  }
  else Alert << "hypothesis K0_1430_K0_1430Hyp not set!!!" <<endmsg;

}

Hyp1Lh::Hyp1Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsPsi2STo2K2PiGamLh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates())
  ,_K1_1270Hyp(true)
  ,_K1_1400Hyp(false)
  ,_K0_1430_K0_1430Hyp(true)
{
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K1_1270Hyp");

  if (iter !=hypMap.end()){
    _K1_1270Hyp= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1270Hyp <<endmsg;
  }
  else Alert << "hypothesis K1_1270Hyp not set!!!" <<endmsg;

  iter= hypMap.find("K1_1400Hyp");
  if (iter !=hypMap.end()){
    _K1_1400Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1400Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }

  else Alert << "hypothesis K1_1400Hyp not set!!!" <<endmsg;
  iter= hypMap.find("K0_1430_K0_1430Hyp");
  if (iter !=hypMap.end()){
    _K0_1430_K0_1430Hyp= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1270Hyp <<endmsg;
  }
  else Alert << "hypothesis K0_1430_K0_1430Hyp not set!!!" <<endmsg;
}

Hyp1Lh::~Hyp1Lh()
{;
}


complex<double> Hyp1Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K892=theParamVal.ChiTo2K892;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_2_1430=theParamVal.ChiTo2K_2_1430;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTo2K_0_1430=theParamVal.ChiTo2K_0_1430;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400ToK892pi=theParamVal.ChiToK1400ToK892pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K1400ToK892Pi=theParamVal.K1400ToK892Pi;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK892Pi=theParamVal.K_1_1270ToK892Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > K_1_1270ToK_0_1430Pi=theParamVal.K_1_1270ToK_0_1430Pi;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1710_k=theParamVal.ChiTof980_pif1710_k;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1710_pi=theParamVal.ChiTof980_kf1710_pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f980=theParamVal.ChiTof980f980;

  double K892Mass=theParamVal.BwK892.first;
  double K892Width=theParamVal.BwK892.second;
  double K_2_1430Mass=theParamVal.BwK_2_1430.first;
  double K_2_1430Width=theParamVal.BwK_2_1430.second;
  double K_0_1430Mass=theParamVal.BwK_0_1430.first;
  double K_0_1430Width=theParamVal.BwK_0_1430.second;
  double K1400Mass=theParamVal.BwK1400.first;
  double K1400Width=theParamVal.BwK1400.second;
  double K_1_1270Mass=theParamVal.BwK_1_1270.first;
  double K_1_1270Width=theParamVal.BwK_1_1270.second;
  double f980_Mass=theParamVal.Flatf980;
  double f980_gPiPi=theParamVal.Flatf980gPiPi;
  double f980_gKK=theParamVal.Flatf980gKK;

  double f1710_kMass=theParamVal.Bwf1710_k.first;
  double f1710_kWidth=theParamVal.Bwf1710_k.second;



  //Chi_c0 decay to K*(892) K*(892)
  complex<double> ChiTo2K892Amp=chiTo2K892Amp(theData, ChiTo2K892, K892Mass, K892Width);

  //Chi_c0 decay to K2*(1430) K2*(1430)
  complex <double> ChiTo2K_1430Amp(0.,0.);
  if (_K0_1430_K0_1430Hyp) ChiTo2K_1430Amp=chiTo2K_0_Amp(theData, ChiTo2K_0_1430, K_0_1430Mass, K_0_1430Width, K_0_1430Mass, K_0_1430Width);
  else ChiTo2K_1430Amp=chiTo2K_2_Amp(theData, ChiTo2K_2_1430, K_2_1430Mass, K_2_1430Width);


  //Chi_c0 decay to K1*(1400) -> K1*(892) pi0 -> (K pi0) pi0
  complex<double> ChiToK_1_1400Amp(0.,0.);

  if (_K1_1270Hyp){
    ChiToK_1_1400Amp+=chiToK1ToK1piAmp(theData, K_1_1270ToK892Pi, K_1_1270Mass, K_1_1270Width, K892Mass, K892Width);
    ChiToK_1_1400Amp+=chiToK1ToK0piAmp(theData, K_1_1270ToK_0_1430Pi, K_1_1270Mass, K_1_1270Width, K_0_1430Mass, K_0_1430Width);
  }
  if (_K1_1400Hyp){
    ChiToK_1_1400Amp+=chiToK1ToK1piAmp(theData, K1400ToK892Pi, K1400Mass, K1400Width, K892Mass, K892Width);
  }
  
    //Chi_c0 decay to f0(980) f0(1710) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980_pif1710_kAmp=chiTof980_pif0_kAmp(theData, ChiTof980_pif1710_k, f980_Mass, f980_gPiPi,  f980_gKK, f1710_kMass, f1710_kWidth);

  //Chi_c0 decay to f0(980) f0(1710) -> (K K) (pi0 pi0) 
  complex<double> ChiTof980_kf1710_piAmp=chiTof980_kf0_piAmp(theData, ChiTof980_kf1710_pi, f980_Mass, f980_gKK,  f980_gPiPi, f1710_kMass, f1710_kWidth);

  //Chi_c0 decay to f0(980) f0(980) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980f980Amp=chiTof980f980Amp(theData, ChiTof980f980, f980_Mass, f980_gPiPi,  f980_gKK);

  complex<double> result=ChiTo2K892Amp+ChiTo2K_1430Amp+ChiToK_1_1400Amp+ChiTof980_pif1710_kAmp+ChiTof980_kf1710_piAmp+ChiTof980f980Amp;

  return result; 

}



void Hyp1Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  AbsPsi2STo2K2PiGamLh::setMnUsrParams(upar, startVal, errVal);

  setMnUsrParamsDec(upar, startVal, errVal,"ChiGam");
  setMnUsrParamsDec(upar, startVal, errVal,"K892K892");

  if (_K0_1430_K0_1430Hyp) setMnUsrParamsDec(upar, startVal, errVal,"K_0_1430K_0_1430");
  else setMnUsrParamsDec(upar, startVal, errVal,"K_2_1430K_2_1430");

  if (_K1_1270Hyp){ 
    setMnUsrParamsDec(upar, startVal, errVal,"K_1_1270ToK892Pi");
    setMnUsrParamsDec(upar, startVal, errVal,"K_1_1270ToK_0_1430Pi");
  }

  if (_K1_1400Hyp){ 
    setMnUsrParamsDec(upar, startVal, errVal,"K_1_1400K");
  }

  setMnUsrParamsDec(upar, startVal, errVal,"f980_pif1710_k");
  setMnUsrParamsDec(upar, startVal, errVal,"f980_kf1710_pi");
  setMnUsrParamsDec(upar, startVal, errVal,"f980f980");

  //K*(892) BW
  setMnUsrParamsMass(upar, startVal, errVal, "K892");

  //K2*(1430) BW
  if (_K0_1430_K0_1430Hyp) setMnUsrParamsMass(upar, startVal, errVal, "K_0_1430"); 
  else setMnUsrParamsMass(upar, startVal, errVal, "K_2_1430");

  //K1(1400) BW
  if (_K1_1270Hyp){ 
    setMnUsrParamsMass(upar, startVal, errVal, "K_1_1270");
    if (!_K0_1430_K0_1430Hyp) setMnUsrParamsMass(upar, startVal, errVal, "K_0_1430"); 
  }

  if (_K1_1400Hyp){ 
    setMnUsrParamsMass(upar, startVal, errVal, "K_1_1400");
  }


//   setMnUsrParamsMass(upar, startVal, errVal, "f980_pi");
  setMnUsrParamsFlatteMass(upar, startVal, errVal, "f980");
  setMnUsrParamsMass(upar, startVal, errVal, "f1710_k");
}



int Hyp1Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_1430States;
  if (_K0_1430_K0_1430Hyp) JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();
  else JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();

  size_t K1400ToKstPiStatesSize=0;
  size_t addResSize=0;
  if (_K1_1270Hyp){
    K1400ToKstPiStatesSize+=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates().size();
    K1400ToKstPiStatesSize+=_Psi2STo2K2PiGamStatesPtr->K1ToK0PiStates().size();
    if(!_K0_1430_K0_1430Hyp) addResSize=1;
  }
  if (_K1_1400Hyp){
    K1400ToKstPiStatesSize+=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates().size();
  }

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();

  if (par.size() < JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2
      +K1400ToKstPiStatesSize*2+JPCLSChiTof0f0States.size()*6+13+addResSize) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be more than " << JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2+K1400ToKstPiStatesSize*2+JPCLSChiTof0f0States.size()*6+13+addResSize << endmsg;
    exit(1);
  }



  int counter=AbsPsi2STo2K2PiGamLh::setFitParamVal(theParamVal, par);

  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "ChiGam");

  //K*(982) K*(982)    amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "K892K892");

  //K2*(1430) K2*(1430)    amplitude params
  if (_K0_1430_K0_1430Hyp) counter=setFitParamValDec(theParamVal, par, counter, "K_0_1430K_0_1430");
  else counter=setFitParamValDec(theParamVal, par, counter, "K_2_1430K_2_1430");

  //K1(1400) K    amplitude params
  if (_K1_1270Hyp){
    counter=setFitParamValDec(theParamVal, par, counter, "K_1_1270ToK892Pi");
    counter=setFitParamValDec(theParamVal, par, counter, "K_1_1270ToK_0_1430Pi");
  }

  if (_K1_1400Hyp){
    counter=setFitParamValDec(theParamVal, par, counter, "K_1_1400K");
  }

  //f980 f1710    amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "f980_pif1710_k");
  counter=setFitParamValDec(theParamVal, par, counter, "f980_kf1710_pi");

  //f980 f980    amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "f980f980");

  // K*892 mass 
  counter=setFitParamValMass(theParamVal, par, counter, "K892");

  //K2*(1430) mass
  if (_K0_1430_K0_1430Hyp)  counter=setFitParamValMass(theParamVal, par, counter, "K_0_1430");
  else counter=setFitParamValMass(theParamVal, par, counter, "K_2_1430");

  //K1(1400) mass
  if (_K1_1270Hyp){
    counter=setFitParamValMass(theParamVal, par, counter, "K_1_1270");
    if (!_K0_1430_K0_1430Hyp) counter=setFitParamValMass(theParamVal, par, counter, "K_0_1430");
  }
  if (_K1_1400Hyp){
    counter=setFitParamValMass(theParamVal, par, counter, "K_1_1400");
  }


  //f980 mass
//   counter=setFitParamValMass(theParamVal, par, counter, "f980_pi");
  counter=setFitParamFlatteMass(theParamVal, par, counter, "f980");

  counter=setFitParamValMass(theParamVal, par, counter, "f1710_k");

  return counter;
}

void Hyp1Lh::print(std::ostream& os) const{
  os << "Hyp1Lh::print\n";
}

void Hyp1Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  AbsPsi2STo2K2PiGamLh::printCurrentFitResult(theParamVal);

//  print fit paramss
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_2_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_0_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1270ToK0PiStates=_Psi2STo2K2PiGamStatesPtr->K1ToK0PiStates();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();

  for ( itJPCLS=JPCLSChiTo2K892States.begin(); itJPCLS!=JPCLSChiTo2K892States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "K892K892" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTo2K892[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  if (_K0_1430_K0_1430Hyp){
    for ( itJPCLS=JPCLSChiTo2K_0_1430States.begin(); itJPCLS!=JPCLSChiTo2K_0_1430States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_0_1430K_0_1430" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.ChiTo2K_0_1430[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }
  }
  else {
    for ( itJPCLS=JPCLSChiTo2K_2_1430States.begin(); itJPCLS!=JPCLSChiTo2K_2_1430States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_2_1430K_2_1430" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.ChiTo2K_2_1430[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }
  }
  if (_K1_1270Hyp){
    for ( itJPCLS=JPCLSK1400ToKst1PiStates.begin(); itJPCLS!=JPCLSK1400ToKst1PiStates.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_1_1270ToK892Pi" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.K_1_1270ToK892Pi[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }
    for ( itJPCLS=JPCLSK1270ToK0PiStates.begin(); itJPCLS!=JPCLSK1270ToK0PiStates.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_1_1270ToK_0_1430Pi" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.K_1_1270ToK_0_1430Pi[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }

  }
  if (_K1_1400Hyp){
    for ( itJPCLS=JPCLSK1400ToKst1PiStates.begin(); itJPCLS!=JPCLSK1400ToKst1PiStates.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "K_1_1400K" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.K1400ToK892Pi[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }
  }
  
  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "f980_pif1710_k" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_pif1710_k[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "f980_kf1710_pi" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_kf1710_pi[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "f980f980" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980f980[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  DebugMsg<< "K892:" << endmsg;
  std::pair<double, double> tmpParamK892=theParamVal.BwK892;
  DebugMsg <<"\t mass:" << tmpParamK892.first <<"\t width:" << tmpParamK892.second  << endmsg;
  
  if (_K0_1430_K0_1430Hyp){
    DebugMsg<< "K_0_1430:" << endmsg;
    std::pair<double, double> tmpParamK_0_1430=theParamVal.BwK_0_1430;
    DebugMsg <<"\t mass:" << tmpParamK_0_1430.first <<"\t width:" << tmpParamK_0_1430.second  << endmsg;
  }
  else{
    DebugMsg<< "K_2_1430:" << endmsg;
    std::pair<double, double> tmpParamK_2_1430=theParamVal.BwK_2_1430;
    DebugMsg <<"\t mass:" << tmpParamK_2_1430.first <<"\t width:" << tmpParamK_2_1430.second  << endmsg;
  }

  if ( _K1_1270Hyp){
    DebugMsg<< "K_1_1270:" << endmsg;
    std::pair<double, double> tmpParamK_1_1270=theParamVal.BwK_1_1270;
    DebugMsg <<"\t mass:" << tmpParamK_1_1270.first <<"\t width:" << tmpParamK_1_1270.second  << endmsg;
    if(!_K0_1430_K0_1430Hyp){
      DebugMsg<< "K_0_1430:" << endmsg;
      std::pair<double, double> tmpParamK_0_1430=theParamVal.BwK_0_1430;
      DebugMsg <<"\t mass:" << tmpParamK_0_1430.first <<"\t width:" << tmpParamK_0_1430.second  << endmsg;
    }
  }

  if ( _K1_1400Hyp){
    DebugMsg<< "K_1_1400:" << endmsg;
    std::pair<double, double> tmpParamK_1_1400=theParamVal.BwK1400;
    DebugMsg <<"\t mass:" << tmpParamK_1_1400.first <<"\t width:" << tmpParamK_1_1400.second  << endmsg;
  }
//     DebugMsg<< "f980_pi:" << endmsg;
//     std::pair<double, double> tmpParamf980_pi=theParamVal.Bwf980_pi;
//     DebugMsg <<"\t mass:" << tmpParamf980_pi.first <<"\t width:" << tmpParamf980_pi.second  << endmsg;

    DebugMsg<< "f980:" << endmsg;
    double tmpParamf980Mass=theParamVal.Flatf980;
    double tmpParamf980gPiPi=theParamVal.Flatf980gPiPi;
    double tmpParamf980gKK=theParamVal.Flatf980gKK;
    DebugMsg <<"\t mass:" << tmpParamf980Mass <<"\t gPiPi:" << tmpParamf980gPiPi  <<"\t gKK:" << tmpParamf980gKK << endmsg;


    DebugMsg<< "f1710_k:" << endmsg;
    std::pair<double, double> tmpParamf1710_k=theParamVal.Bwf1710_k;
    DebugMsg <<"\t mass:" << tmpParamf1710_k.first <<"\t widthKK:" << tmpParamf1710_k.second << endmsg;



}

void Hyp1Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  AbsPsi2STo2K2PiGamLh::dumpCurrentResult(os, theParamVal, suffix);

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_2_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_0_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1270ToK0PiStates=_Psi2STo2K2PiGamStatesPtr->K1ToK0PiStates();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();

  std::string tmpStringDec("0");

  for ( itJPCLS=JPCLSChiTo2K892States.begin(); itJPCLS!=JPCLSChiTo2K892States.end(); ++itJPCLS){
    tmpStringDec=(*itJPCLS)->name()+"K892K892"+suffix;
    std::pair<double, double> tmpParam=theParamVal.ChiTo2K892[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }
  
  if (_K0_1430_K0_1430Hyp){
    for ( itJPCLS=JPCLSChiTo2K_0_1430States.begin(); itJPCLS!=JPCLSChiTo2K_0_1430States.end(); ++itJPCLS){
      tmpStringDec=(*itJPCLS)->name()+"K_0_1430K_0_1430"+suffix;
      std::pair<double, double> tmpParam=theParamVal.ChiTo2K_0_1430[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
  }
  else {
    for ( itJPCLS=JPCLSChiTo2K_2_1430States.begin(); itJPCLS!=JPCLSChiTo2K_2_1430States.end(); ++itJPCLS){
      tmpStringDec=(*itJPCLS)->name()+"K_2_1430K_2_1430"+suffix;
      std::pair<double, double> tmpParam=theParamVal.ChiTo2K_2_1430[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
  }

  if (_K1_1270Hyp){

    for ( itJPCLS=JPCLSK1400ToKst1PiStates.begin(); itJPCLS!=JPCLSK1400ToKst1PiStates.end(); ++itJPCLS){
      tmpStringDec=(*itJPCLS)->name()+"K_1_1270ToK892Pi"+suffix;
      std::pair<double, double> tmpParam=theParamVal.K_1_1270ToK892Pi[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }

    for ( itJPCLS=JPCLSK1270ToK0PiStates.begin(); itJPCLS!=JPCLSK1270ToK0PiStates.end(); ++itJPCLS){
      tmpStringDec=(*itJPCLS)->name()+"K_1_1270ToK_0_1430Pi"+suffix;
      std::pair<double, double> tmpParam=theParamVal.K_1_1270ToK_0_1430Pi[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;  
    }

  }

  if (_K1_1400Hyp){
    for ( itJPCLS=JPCLSK1400ToKst1PiStates.begin(); itJPCLS!=JPCLSK1400ToKst1PiStates.end(); ++itJPCLS){
      tmpStringDec=(*itJPCLS)->name()+"K_1_1400K"+suffix;
      std::pair<double, double> tmpParam=theParamVal.K1400ToK892Pi[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
  }
  
  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    tmpStringDec=(*itJPCLS)->name()+"f980_pif1710_k"+suffix;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_pif1710_k[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    tmpStringDec=(*itJPCLS)->name()+"f980_kf1710_pi"+suffix;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980_kf1710_pi[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

  for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
    tmpStringDec=(*itJPCLS)->name()+"f980f980"+suffix;
    std::pair<double, double> tmpParam=theParamVal.ChiTof980f980[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

  std::string tmpStringRes="K892mass"+suffix;
  std::pair<double, double> tmpParamRes=theParamVal.BwK892;
  os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;

  if (_K0_1430_K0_1430Hyp){
    tmpStringRes="K_0_1430mass"+suffix;
    tmpParamRes=theParamVal.BwK_0_1430;
    os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;
  }
  else{
    tmpStringRes="K_2_1430mass"+suffix;
    tmpParamRes=theParamVal.BwK_2_1430;
    os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;
  }

  if (_K1_1270Hyp){
    tmpStringRes="K_1_1270mass"+suffix;
    tmpParamRes=theParamVal.BwK_1_1270;
    os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;

    if (!_K0_1430_K0_1430Hyp){
      tmpStringRes="K_0_1430mass"+suffix;
      tmpParamRes=theParamVal.BwK_0_1430;
      os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;
    }
  }

  if (_K1_1400Hyp){
    tmpStringRes="K_1_1400mass"+suffix;
    tmpParamRes=theParamVal.BwK1400;
    os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;
  }

  std::string tmpStringFlatMass="FlatMf980"+suffix;
  std::string tmpStringFlatPiPi="Flatf980gPiPi"+suffix;
  std::string tmpStringFlatKK="Flatf980gKK"+suffix;
  double tmpParamFlatMass = theParamVal.Flatf980;
  double tmpParamFlatPiPi=theParamVal.Flatf980gPiPi;
  double tmpParamFlatKK=theParamVal.Flatf980gKK;
  os << tmpStringFlatMass << "\t" << tmpParamFlatMass <<"\t" << 0.  << std::endl;
  os << tmpStringFlatPiPi << "\t" << tmpParamFlatPiPi <<"\t" << 0.  << std::endl;
  os << tmpStringFlatKK << "\t" << tmpParamFlatKK <<"\t" << 0.  << std::endl;

  tmpStringRes="f1710_kmass"+suffix;
  tmpParamRes=theParamVal.Bwf1710_k;
  os << tmpStringRes << "\t" << tmpParamRes.first <<"\t" << tmpParamRes.second  << std::endl;

}
