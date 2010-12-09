#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

Hyp3Lh::Hyp3Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates, const std::map<const std::string, bool>& hypMap ) :
  Hyp2Lh(theEvtList,theStates, hypMap)
  ,_sigmaf1710Hyp(true)
  ,_sigmaf980Hyp(true)
  ,_disableHyp3(false)
{
 
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("sigmaf1710Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf1710Hyp= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf1710Hyp <<endmsg;
  }
  else Alert << "hypothesis sigmaf1710Hyp3 not set!!!" <<endmsg;

 iter= hypMap.find("sigmaf980Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf980Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf980Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis sigmaf980Hyp not set!!!" <<endmsg;

  if (!_sigmaf980Hyp && !_sigmaf1710Hyp) _disableHyp3=true; 
}

Hyp3Lh::Hyp3Lh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  Hyp2Lh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates(), hypMap)
  ,_sigmaf1710Hyp(true)
  ,_sigmaf980Hyp(true)
  ,_disableHyp3(false)
{

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("sigmaf1710Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf1710Hyp= iter->second;
   _hypMap[iter->first]= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf1710Hyp <<endmsg;
  }
  else Alert << "hypothesis sigmaf1710Hyp3 not set!!!" <<endmsg;

 iter= hypMap.find("sigmaf980Hyp3");

  if (iter !=hypMap.end()){
    _sigmaf980Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _sigmaf980Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis sigmaf980Hyp not set!!!" <<endmsg;

  if (!_sigmaf980Hyp && !_sigmaf1710Hyp) _disableHyp3=true; 
}

Hyp3Lh::~Hyp3Lh()
{;
}


complex<double> Hyp3Lh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  complex<double> result=Hyp2Lh::chi0DecAmps(theParamVal, theData);

  if(_disableHyp3) return result;
 
  double sigmaMass=theParamVal.BwSigma.first;
  double sigmaWidth=theParamVal.BwSigma.second;

  if (_sigmaf980Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf980=theParamVal.ChiToSigmaf980;
 
    double f980_Mass=theParamVal.Flatf980;
    double f980_gPiPi=theParamVal.Flatf980gPiPi;
    double f980_gKK=theParamVal.Flatf980gKK;
    result+= chiTof980_kf0_piAmp(theData, ChiToSigmaf980, f980_Mass, f980_gKK,  f980_gPiPi, sigmaMass, sigmaWidth);
  }

  if (_sigmaf1710Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToSigmaf1710=theParamVal.ChiToSigmaf1710;
    double f1710_kMass=theParamVal.Bwf1710_k.first;
    double f1710_kWidth=theParamVal.Bwf1710_k.second;
    result+=chiTof0_pif0_kAmp(theData, ChiToSigmaf1710, sigmaMass, sigmaWidth,  f1710_kMass, f1710_kWidth);
  }  


  return result;
}


void Hyp3Lh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  Hyp2Lh::setMnUsrParams(upar, startVal, errVal);
  if(_disableHyp3) return;


  if(_sigmaf980Hyp) setMnUsrParamsDec(upar, startVal, errVal,"sigmaf980");
  if(_sigmaf1710Hyp) setMnUsrParamsDec(upar, startVal, errVal,"sigmaf1710");

  setMnUsrParamsMass(upar, startVal, errVal, "sigma");   
}



int Hyp3Lh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{

//   std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();

//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
//   if(_K0_1430_K0_1430Hyp) JPCLSChiTo2K_1430States=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();

//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
//   std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
  
//   if (par.size()< JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2
//       +JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*10+17) {
//     Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
//       "\t it should be more than" << JPCLSPsiToChiGam.size()*2-1+JPCLSChiTo2K892States.size()*2+JPCLSChiTo2K_1430States.size()*2+JPCLSK1400ToKst1PiStates.size()*2+JPCLSChiTof0f0States.size()*10+17 << endmsg;
//     exit(1);
//   } 
  

  int counter=Hyp2Lh::setFitParamVal(theParamVal, par);
  if(_disableHyp3) return counter;


  //sigma f1710    amplitude params
  if(_sigmaf980Hyp) counter=setFitParamValDec(theParamVal, par, counter, "sigmaf980");
  if(_sigmaf1710Hyp) counter=setFitParamValDec(theParamVal, par, counter, "sigmaf1710");
  
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

    if(_sigmaf980Hyp){
      for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
	DebugMsg<< (*itJPCLS)->name()<< "sigmaf980" << endmsg;
	std::pair<double, double> tmpParam=theParamVal.ChiToSigmaf980[(*itJPCLS)];
	DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
	
      }
    }

    if(_sigmaf1710Hyp){
      for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
	DebugMsg<< (*itJPCLS)->name()<< "sigmaf1710" << endmsg;
	std::pair<double, double> tmpParam=theParamVal.ChiToSigmaf1710[(*itJPCLS)];
	DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
	
      }
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

    if(_sigmaf980Hyp){    
      for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
	std::string tmpStringDec=(*itJPCLS)->name()+"sigmaf980"+suffix;
	
	std::pair<double, double> tmpParam=theParamVal.ChiToSigmaf980[(*itJPCLS)];
	os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
      }
    }

    if(_sigmaf1710Hyp){    
      for ( itJPCLS=JPCLSChiTof0f0States.begin(); itJPCLS!=JPCLSChiTof0f0States.end(); ++itJPCLS){
	std::string tmpStringDec=(*itJPCLS)->name()+"sigmaf1710"+suffix;
	
	std::pair<double, double> tmpParam=theParamVal.ChiToSigmaf1710[(*itJPCLS)];
	os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
      }
    }
    
    std::string tmpStringRes="sigmamass"+suffix;
    std::pair<double, double> tmpParamSigma=theParamVal.BwSigma;
    os << tmpStringRes << "\t" << tmpParamSigma.first <<"\t" << tmpParamSigma.second  << std::endl;
  }
}
