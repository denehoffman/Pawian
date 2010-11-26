#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp3Lh::Hyp3Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp2Lh(theEvtList,theStates, hypMap)
  ,_disableHyp3(false)
{
 
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("disableHyp3");

  if (iter !=hypMap.end()){
    _disableHyp3= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _disableHyp3 <<endmsg;
  }
  else Alert << "hypothesis " << iter->first << "not set!!!" <<endmsg;

}

Hyp3Lh::Hyp3Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp2Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
  ,_disableHyp3(false)
{
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("disableHyp3");

  if (iter !=hypMap.end()){
    _disableHyp3= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _disableHyp3 <<endmsg;
  }
  else Alert << "hypothesis " << iter->first << "not set!!!" <<endmsg;
}

Hyp3Lh::~Hyp3Lh()
{;
}


complex<double> Hyp3Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp2Lh::chi0DecAmps(theParamVal, theData);
  if(_disableHyp3) return result;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710=theParamVal.ChiToSigmaf1710;
  double sigmaMass=theParamVal.BwSigma.first;
  double sigmaWidth=theParamVal.BwSigma.second;
  double f1710_kMass=theParamVal.Bwf1710_k.first;
  double f1710_kWidth=theParamVal.Bwf1710_k.second;

  //Chi_c0 decay to sigma f0(1710) -> (pi0 pi0) (K K)
  complex<double> ChiToSigma_pif1710_kAmp=chiTof0_pif0_kAmp(theData, ChiToSigmaf1710, sigmaMass, sigmaWidth,  f1710_kMass, f1710_kWidth);
  
  result+=ChiToSigma_pif1710_kAmp;

  return result;
}


void Hyp3Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  Hyp2Lh::setMnUsrParams(upar, startVal, errVal);

  if(!_disableHyp3){
    setMnUsrParamsDec(upar, startVal, errVal,"sigmaf1710");
    setMnUsrParamsMass(upar, startVal, errVal, "sigma"); 
  } 
}



int Hyp3Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
  if(_K0_1430_K0_1430Hyp) JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  
  if (par.size()< JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2
      +JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*10+17) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be more than" << JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2+JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*10+17 << endmsg;
    exit(1);
  }  

  int counter=Hyp2Lh::setFitParamVal(theParamVal, par);
  if(_disableHyp3) return counter;

  //sigma f1710    amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "sigmaf1710");

  counter=setFitParamValMass(theParamVal, par, counter, "sigma");  

  return counter;
}

void Hyp3Lh::print(std::ostream& os) const{
  os << "Hyp3Lh::print\n";
}

void Hyp3Lh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  Hyp2Lh::printCurrentFitResult(theParamVal);

  if(!_disableHyp3){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< "sigmaf1710" << endmsg;
      std::pair<double, double> tmpParam=theParamVal.ChiToSigmaf1710[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
      
    }
    
    DebugMsg<< "sigma:" << endmsg;
    std::pair<double, double> tmpParamSigma=theParamVal.BwSigma;
    DebugMsg <<"\t mass:" << tmpParamSigma.first <<"\t width:" << tmpParamSigma.second  << endmsg;
  }
}

void Hyp3Lh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{
 
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }
  Hyp2Lh::dumpCurrentResult(os, theParamVal, suffix);

  if(!_disableHyp3){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
    
    for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
      std::string tmpStringDec=(*itJPCLS)->name()+"sigmaf1710"+suffix;
      
      std::pair<double, double> tmpParam=theParamVal.ChiToSigmaf1710[(*itJPCLS)];
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }
    
    std::string tmpStringRes="sigmamass"+suffix;
    std::pair<double, double> tmpParamSigma=theParamVal.BwSigma;
    os << tmpStringRes << "\t" << tmpParamSigma.first <<"\t" << tmpParamSigma.second  << std::endl;
  }
}
