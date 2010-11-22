#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp4Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp4Lh::Hyp4Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, bool K1_1270Hyp, bool K0_1430_K0_1430Hyp, bool disableHyp3) :
  Hyp3Lh(theEvtList,theStates, K1_1270Hyp, K0_1430_K0_1430Hyp, disableHyp3)  
{
}

Hyp4Lh::Hyp4Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, bool K1_1270Hyp, bool K0_1430_K0_1430Hyp, bool disableHyp3 ) :
  Hyp3Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), K1_1270Hyp, K0_1430_K0_1430Hyp, disableHyp3)
{
}

Hyp4Lh::~Hyp4Lh()
{;
}


complex<double> Hyp4Lh::calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;
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
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980f2200=theParamVal.ChiTof980f2200;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710=theParamVal.ChiToSigmaf1710;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_pif1370_k=theParamVal.ChiTof980_pif1370_k;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiTof980_kf1370_pi=theParamVal.ChiTof980_kf1370_pi;

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

  double f2200Mass=theParamVal.Bwf2200.first;
  double f2200Width=theParamVal.Bwf2200.second;

  double sigmaMass=theParamVal.BwSigma.first;
  double sigmaWidth=theParamVal.BwSigma.second;

  double f1370_kMass=theParamVal.Bwf1370_k.first;
  double f1370_kWidth=theParamVal.Bwf1370_k.second;
  double f1710PiPiWidth=theParamVal.f1710PiPiWidth;

  complex<double> PsiDecAmpTmpls(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
  for ( itPsi=PsiToChiGam.begin(); itPsi!=PsiToChiGam.end(); ++itPsi){
    
    boost::shared_ptr<const JPCLS> PsiState=itPsi->first;

    double thePsiMag=itPsi->second.first;
    double thePsiPhi=itPsi->second.second;
    complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));


    PsiDecAmpTmpls+=thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
      *Clebsch(PsiState->L,0,PsiState->S, -lamGam, PsiState->J, -lamGam)
      *Clebsch(0, 0, 1, -lamGam, PsiState->S, -lamGam);
  }

  //Chi_c0 decay to K*(892) K*(892)
  complex<double> ChiTo2K892Amp=chiTo2K892Amp(theData, ChiTo2K892, K892Mass, K892Width);

  //Chi_c0 decay to K2*(1430) K2*(1430)
  complex <double> ChiTo2K_1430Amp(0.,0.);
  if (_K0_1430_K0_1430Hyp) ChiTo2K_1430Amp=chiTo2K_0_Amp(theData, ChiTo2K_0_1430, K_0_1430Mass, K_0_1430Width);
  else ChiTo2K_1430Amp=chiTo2K_2_Amp(theData, ChiTo2K_2_1430, K_2_1430Mass, K_2_1430Width);

  //Chi_c0 decay to K1*(1400) -> K1*(892) pi0 -> (K pi0) pi0
  complex<double> ChiToK_1_1400Amp(0.,0.);

  if (_K1_1270Hyp){
    ChiToK_1_1400Amp=chiToK1ToK1piAmp(theData, K_1_1270ToK892Pi, K_1_1270Mass, K_1_1270Mass, K892Mass, K892Width);
    ChiToK_1_1400Amp+=chiToK1ToK0piAmp(theData, K_1_1270ToK_0_1430Pi, K_1_1270Mass, K_1_1270Width, K_0_1430Mass, K_0_1430Width);
  }
  else ChiToK_1_1400Amp=chiToK1ToK1piAmp(theData, K1400ToK892Pi, K1400Mass, K1400Width, K892Mass, K892Width);


  //Chi_c0 decay to f0(980) f0(1710) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980_pif1710_kAmp=chiTof980_pif0_kAmp(theData, ChiTof980_pif1710_k, f980_Mass, f980_gPiPi,  f980_gKK, f1710_kMass, f1710_kWidth);
  complex<double> ChiTof980_kf1710_piAmp=chiTof980_kf0_piAmp(theData, ChiTof980_kf1710_pi, f980_Mass, f980_gKK,  f980_gPiPi, f1710_kMass, f1710PiPiWidth);

  //Chi_c0 decay to f0(980) f0(980) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980f980Amp=chiTof980f980Amp(theData, ChiTof980f980, f980_Mass, f980_gPiPi,  f980_gKK);

  //Chi_c0 decay to f0(980) f0(2200) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980f2200Amp=chiTof980_pif0_kAmp(theData, ChiTof980f2200, f980_Mass, f980_gPiPi,  f980_gKK, f2200Mass, f2200Width);

  //Chi_c0 decay to sigma f0(1710) -> (pi0 pi0) (K K)
  complex<double> ChiToSigma_pif1710_kAmp(0.,0.); 
  if(!_disableHyp3){
    ChiToSigma_pif1710_kAmp=chiTof0_pif0_kAmp(theData, ChiToSigmaf1710, sigmaMass, sigmaWidth,  f1710_kMass, f1710_kWidth);
  }

  //Chi_c0 decay to f0(980) f0(1370) -> (pi0 pi0) (K K) 
  complex<double> ChiTof980_pif1370_kAmp=chiTof980_pif0_kAmp(theData, ChiTof980_pif1370_k, f980_Mass, f980_gPiPi,  f980_gKK, f1370_kMass, f1370_kWidth);
  complex<double> ChiTof980_kf1370_piAmp=chiTof980_kf0_piAmp(theData, ChiTof980_kf1370_pi, f980_Mass, f980_gKK,  f980_gPiPi, f1370_kMass, f1370_kWidth);

  result+=conj(theData->DfPsi[1][Minit][-lamGam])*PsiDecAmpTmpls
    *(ChiTo2K892Amp+ChiTo2K_1430Amp+ChiToK_1_1400Amp+ChiTof980_pif1710_kAmp+ChiTof980_kf1710_piAmp+ChiTof980f980Amp+ChiTof980f2200Amp+ChiToSigma_pif1710_kAmp+ChiTof980_pif1370_kAmp+ChiTof980_kf1370_piAmp);



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
