#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsPsi2STo2K2PiGamLh::AbsPsi2STo2K2PiGamLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates) :
  _Psi2STo2K2PiGamEvtListPtr(theEvtList),
  _Psi2STo2K2PiGamStatesPtr(theStates)
{
  _evtDataVec=_Psi2STo2K2PiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2STo2K2PiGamEvtListPtr->getMcVecs();
}

AbsPsi2STo2K2PiGamLh::AbsPsi2STo2K2PiGamLh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theAbsPsi2STo2K2PiGamLhPtr):
  _Psi2STo2K2PiGamEvtListPtr(theAbsPsi2STo2K2PiGamLhPtr->getEventList()),
  _Psi2STo2K2PiGamStatesPtr(theAbsPsi2STo2K2PiGamLhPtr->getPsi2STo2K2PiGamStates())
{
  _evtDataVec=_Psi2STo2K2PiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2STo2K2PiGamEvtListPtr->getMcVecs();
}

AbsPsi2STo2K2PiGamLh::~AbsPsi2STo2K2PiGamLh()
{
}

double AbsPsi2STo2K2PiGamLh::calcLogLh(const Psi2STo2K2PiGamData::fitParamVal& theParamVal){
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData*>::iterator iterm;
  for (iterm=_evtMCVec.begin(); iterm!=_evtMCVec.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log10(LH_mc/_evtMCVec.size());

  logLH=_evtDataVec.size()/2.*(LH_mc/_evtMCVec.size()-1)*(LH_mc/_evtMCVec.size()-1)
    -logLH_data
    +_evtDataVec.size()*logLH_mc_Norm;

  Info << "current LH = " << logLH << endmsg;

 return logLH;

}

double AbsPsi2STo2K2PiGamLh::calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const Psi2STo2K2PiGamData::fitParamVal& theParamVal){
  double phaseSpaceVal=theParamVal.phaseSpace;
  
  complex<double> theDecAmp=chi0DecAmps(theParamVal, theData);

  Spin Psi2SM=1;
  Spin GamM=1;
  complex<double> AmpPsi2SMpGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

  Psi2SM=1;
  GamM=-1;
  complex<double> AmpPsi2SMpGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

  Psi2SM=-1;
  GamM=1; 
  complex<double> AmpPsi2SMmGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

  Psi2SM=-1;
  GamM=-1; 
  complex<double> AmpPsi2SMmGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData)*theDecAmp;

//   DebugMsg << "AmpPsi2SMp " << AmpPsi2SMp << endmsg;
  
  double result=norm(AmpPsi2SMpGp)+norm(AmpPsi2SMpGm)+norm(AmpPsi2SMmGp)+norm(AmpPsi2SMmGm)+phaseSpaceVal;

  return result;  
}
  
complex<double> AbsPsi2STo2K2PiGamLh::calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;

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

  complex<double> result=conj(theData->DfPsi[1][Minit][-lamGam])*PsiDecAmpTmpls;
  return result; 
}



void AbsPsi2STo2K2PiGamLh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  double phaseSpaceVal=startVal.phaseSpace;
  double phaseSpaceErr=startVal.phaseSpace;
  std::string phaseSpaceStr="PhaseSpace";

  double phaseSpaceMax=phaseSpaceVal+3*phaseSpaceErr;
  if ( phaseSpaceMax > 1. ) phaseSpaceMax=1.;
 
  upar.Add(phaseSpaceStr, phaseSpaceVal, phaseSpaceErr, 0., phaseSpaceMax);
}

int AbsPsi2STo2K2PiGamLh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
 int counter=0;
 theParamVal.phaseSpace=par[counter];
 counter++;
 return counter;
}

int AbsPsi2STo2K2PiGamLh::setFitParamValDec(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const{

  int resultCount=counter;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > currentStates;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >* currentMap=0;


  if (key=="ChiGam"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
    currentMap= &theParamVal.PsiToChiGam;
  }
  else if (key=="K892K892"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
    currentMap= &theParamVal.ChiTo2K892;
  }
  else if (key=="K_0_1430K_0_1430"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_0_States();
    currentMap= &theParamVal.ChiTo2K_0_1430;
  }
  else if (key=="K_2_1430K_2_1430"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTo2K_2_1430States();
    currentMap= &theParamVal.ChiTo2K_2_1430;
  }
  else if (key=="K_0_1430K_2_1430"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiToK0K2_States();
    currentMap= &theParamVal.ChiToK_0_1430_K_2_1430;
  }
  else if (key=="K_1_1400K"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
    currentMap= &theParamVal.K1400ToK892Pi;
  }
  else if (key=="K_1_1270ToK892Pi"){
    currentStates=_Psi2STo2K2PiGamStatesPtr->K1400ToKst1PiStates();
    currentMap= &theParamVal.K_1_1270ToK892Pi;
  }
  else if (key=="K_1_1270ToK_0_1430Pi"){
    currentStates=_Psi2STo2K2PiGamStatesPtr->K1ToK0PiStates();
    currentMap= &theParamVal.K_1_1270ToK_0_1430Pi;
  }
  else if (key=="f980_pif1710_k"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiTof980_pif1710_k;
  }
  else if (key=="f980_kf1710_pi"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiTof980_kf1710_pi;
  }
  else if (key=="f980f980"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiTof980f980;
  }
  else if (key=="f980f2200"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiTof980f2200;
  }
  else if (key=="sigmaf980"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiToSigmaf980;
  }
  else if (key=="sigmaf1710"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiToSigmaf1710;
  }
  else if (key=="sigmaf2200"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiToSigmaf2200;
  }
  else if (key=="f980_pif1370_k"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiTof980_pif1370_k;
  }
  else if (key=="f980_kf1370_pi"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTof0f0States();
    currentMap= &theParamVal.ChiTof980_kf1370_pi;
  }
  else if (key=="K_0_2400KToKf980K"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    currentMap= &theParamVal.ChiToK_0_2400ToKf980;
  }
  else if (key=="K_0_1430K_0_1950"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    currentMap= &theParamVal.ChiToK_0_1430K_0_1950;
  }
  else if (key=="K892K_1_1680"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
    currentMap= &theParamVal.ChiToK892K1680;
  }
  else if (key=="K892K_1_2300"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiTo2K892States();
    currentMap= &theParamVal.ChiToK892K2300;
  }
  else if (key=="K_0_1460ToK892Pi"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    currentMap= &theParamVal.K_0_1460ToK892Pi;
  }
  else if (key=="K_0_1460ToK_0_1430Pi"){ 
    currentStates=_Psi2STo2K2PiGamStatesPtr->ChiToK0K0States();
    currentMap= &theParamVal.K_0_1460ToK_0_1430Pi;
  }

  else { Alert << "Key: " << key << " not supported for setting up the fit amplitude parameters!!!" << endmsg;
    exit(1);
  }

  for ( itJPCLS=currentStates.begin(); itJPCLS!=currentStates.end(); ++itJPCLS){
    double mag=par[resultCount];
    resultCount++;
    double phi=par[resultCount];
    resultCount++;
    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    (*currentMap)[(*itJPCLS)]=tmpParameter;
  }

  return resultCount;
}


int AbsPsi2STo2K2PiGamLh::setFitParamValMass(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const{

  int resultCount=counter;

  double currentMass=par[resultCount];
  resultCount++;
  double currentWidth=par[resultCount];
  resultCount++;

  std::pair <double,double> thePair=make_pair(currentMass, currentWidth);

  if (key=="K892"){
    theParamVal.BwK892=thePair;
  }
  else if (key=="K_0_1430"){
    theParamVal.BwK_0_1430=thePair;
  }
  else if (key=="K_2_1430"){
    theParamVal.BwK_2_1430=thePair;
  }
  else if (key=="K_1_1400"){
    theParamVal.BwK1400=thePair;
  }
  else if (key=="K_1_1270"){
    theParamVal.BwK_1_1270=thePair;
  }
  else if (key=="f980_pi"){
    theParamVal.Bwf980_pi=thePair;
  }
  else if (key=="f1710_k"){
    theParamVal.Bwf1710_k=thePair;
  }
  else if (key=="f2200"){
    theParamVal.Bwf2200=thePair;
  }
  else if (key=="sigma"){
    theParamVal.BwSigma=thePair;
  }
  else if (key=="f1370_k"){
    theParamVal.Bwf1370_k=thePair;
  }
  else if (key=="K_0_2400"){
    theParamVal.BwK_0_2400=thePair;
  }
  else if (key=="K_0_1950"){
    theParamVal.BwK_0_1950=thePair;
  }
  else if (key=="K_1_1680"){
    theParamVal.BwK_1_1680=thePair;
  }
  else if (key=="K_1_2300"){
    theParamVal.BwK_1_2300=thePair;
  }
  else if (key=="K_0_1460"){
    theParamVal.BwK_0_1460=thePair;
  }
  else{
    Alert << "Key: " << key << " not supported for setting up the fit parameter for masses!!!" << endmsg;
    exit(1);
  }

  return resultCount;
}

int AbsPsi2STo2K2PiGamLh::setFitParamFlatteMass(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const{

  int resultCount=counter;

  if (key=="f980"){
    theParamVal.Flatf980=par[resultCount];
    resultCount++;
    theParamVal.Flatf980gPiPi=par[resultCount];
    resultCount++;
    theParamVal.Flatf980gKK=par[resultCount];
    resultCount++;
  }
  else { Alert << "Key: " << key << " not supported for setting up the Flatte mass parameters!!!" << endmsg;
    exit(1);
  }


  return resultCount;
}

void AbsPsi2STo2K2PiGamLh::setMnUsrParamsDec(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal, std::string key){
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > startParams;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > errParams;

  if (key=="ChiGam"){
    startParams=startVal.PsiToChiGam;
    errParams=errVal.PsiToChiGam;
  }
  else if (key=="K892K892"){
    startParams=startVal.ChiTo2K892;
    errParams=errVal.ChiTo2K892;
  }
  else if (key=="K_0_1430K_0_1430"){
    startParams=startVal.ChiTo2K_0_1430;
    errParams=errVal.ChiTo2K_0_1430;
  }
  else if (key=="K_2_1430K_2_1430"){
    startParams=startVal.ChiTo2K_2_1430;
    errParams=errVal.ChiTo2K_2_1430;
  }
  else if (key=="K_0_1430K_2_1430"){
    startParams=startVal.ChiToK_0_1430_K_2_1430;
    errParams=errVal.ChiToK_0_1430_K_2_1430;
  }
  else if (key=="K_1_1400K"){
    startParams=startVal.K1400ToK892Pi;
    errParams=errVal.K1400ToK892Pi;
  }
  else if (key=="K_1_1270ToK892Pi"){
    startParams=startVal.K_1_1270ToK892Pi;
    errParams=errVal.K_1_1270ToK892Pi;
  }
  else if (key=="K_1_1270ToK_0_1430Pi"){
    startParams=startVal.K_1_1270ToK_0_1430Pi;
    errParams=errVal.K_1_1270ToK_0_1430Pi;
  }
  else if (key=="f980_pif1710_k"){
    startParams=startVal.ChiTof980_pif1710_k;
    errParams=errVal.ChiTof980_pif1710_k;
  }
  else if (key=="f980_kf1710_pi"){
    startParams=startVal.ChiTof980_kf1710_pi;
    errParams=errVal.ChiTof980_kf1710_pi;
  }
  else if (key=="f980f980"){
    startParams=startVal.ChiTof980f980;
    errParams=errVal.ChiTof980f980;
  }
  else if (key=="f980f2200"){
    startParams=startVal.ChiTof980f2200;
    errParams=errVal.ChiTof980f2200;
  }
  else if (key=="sigmaf980"){
    startParams=startVal.ChiToSigmaf980;
    errParams=errVal.ChiToSigmaf980;
  }
  else if (key=="sigmaf1710"){
    startParams=startVal.ChiToSigmaf1710;
    errParams=errVal.ChiToSigmaf1710;
  }
  else if (key=="sigmaf2200"){
    startParams=startVal.ChiToSigmaf2200;
    errParams=errVal.ChiToSigmaf2200;
  }
  else if (key=="f980_pif1370_k"){
    startParams=startVal.ChiTof980_pif1370_k;
    errParams=errVal.ChiTof980_pif1370_k;
  }
  else if (key=="f980_kf1370_pi"){
    startParams=startVal.ChiTof980_kf1370_pi;
    errParams=errVal.ChiTof980_kf1370_pi;
  }
  else if (key=="K_0_2400KToKf980K"){
    startParams=startVal.ChiToK_0_2400ToKf980;
    errParams=errVal.ChiToK_0_2400ToKf980;
  }
  else if (key=="K_0_1430K_0_1950"){
    startParams=startVal.ChiToK_0_1430K_0_1950;
    errParams=errVal.ChiToK_0_1430K_0_1950;
  }
  else if (key=="K892K_1_1680"){
    startParams=startVal.ChiToK892K1680;
    errParams=errVal.ChiToK892K1680;
  }
  else if (key=="K892K_1_2300"){
    startParams=startVal.ChiToK892K2300;
    errParams=errVal.ChiToK892K2300;
  }
  else if (key=="K_0_1460ToK892Pi"){
    startParams=startVal.K_0_1460ToK892Pi;
    errParams=errVal.K_0_1460ToK892Pi;
  }
  else if (key=="K_0_1460ToK_0_1430Pi"){
    startParams=startVal.K_0_1460ToK_0_1430Pi;
    errParams=errVal.K_0_1460ToK_0_1430Pi;
  }
  else { Alert << "Key: " << key << " not supported for setting up the MINUIT start amplitude parameters!!!" << endmsg;
    exit(1);
  }

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for (it=startParams.begin(); it!=startParams.end(); ++it){
    boost::shared_ptr<const JPCLS> theJPCLS=it->first;   
    std::pair <double,double> startPair=it->second;
    std::pair <double,double> errPair=errParams.find(theJPCLS)->second;

    //now fill the fitParameterMap
    std::string magStr=theJPCLS->name()+key+"Mag";
    std::string phiStr=theJPCLS->name()+key+"phi";

    double magVal=startPair.first;
    double phiVal=startPair.second;
    
    double magErr=errPair.first;
    double phiErr=errPair.second;
    
    double magMin=magVal-magErr;
    if (magMin<0.) magMin=0.;
    
    upar.Add(magStr, magVal, magErr, magMin, magVal+magErr);
    upar.Add(phiStr, phiVal, phiErr, -3*M_PI, 3*M_PI);

  }

}

void AbsPsi2STo2K2PiGamLh::setMnUsrParamsMass(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal, std::string key){
  pair<double, double> startParams;
  pair<double, double> errParams;

  if (key=="K892"){
    startParams=startVal.BwK892;
    errParams=errVal.BwK892;
  }
  else if (key=="K_0_1430"){
    startParams=startVal.BwK_0_1430;
    errParams=errVal.BwK_0_1430;
  }
  else if (key=="K_2_1430"){
    startParams=startVal.BwK_2_1430;
    errParams=errVal.BwK_2_1430;
  }
  else if (key=="K_1_1400"){
    startParams=startVal.BwK1400;
    errParams=errVal.BwK1400;
  }
  else if (key=="K_1_1270"){
    startParams=startVal.BwK_1_1270;
    errParams=errVal.BwK_1_1270;
  }
  else if (key=="f980_pi"){
    startParams=startVal.Bwf980_pi;
    errParams=errVal.Bwf980_pi;
  }
  else if (key=="f1710_k"){
    startParams=startVal.Bwf1710_k;
    errParams=errVal.Bwf1710_k;
  }
  else if (key=="f2200"){
    startParams=startVal.Bwf2200;
    errParams=errVal.Bwf2200;
  }
  else if (key=="sigma"){
    startParams=startVal.BwSigma;
    errParams=errVal.BwSigma;
  }
  else if (key=="f1370_k"){
    startParams=startVal.Bwf1370_k;
    errParams=errVal.Bwf1370_k;
  }
  else if (key=="K_0_2400"){
    startParams=startVal.BwK_0_2400;
    errParams=errVal.BwK_0_2400;
  }
  else if (key=="K_0_1950"){
    startParams=startVal.BwK_0_1950;
    errParams=errVal.BwK_0_1950;
  }
  else if (key=="K_1_1680"){
    startParams=startVal.BwK_1_1680;
    errParams=errVal.BwK_1_1680;
  }
  else if (key=="K_1_2300"){
    startParams=startVal.BwK_1_2300;
    errParams=errVal.BwK_1_2300;
  }
  else if (key=="K_0_1460"){
    startParams=startVal.BwK_0_1460;
    errParams=errVal.BwK_0_1460;
  }
  else { Alert << "Key: " << key << " not supported for setting up the MINUIT start mass parameters!!!" << endmsg;
    exit(1);
  }

  //now fill the fitParameterMap
  std::string massStr=key+"Mag";
  std::string widthStr=key+"Width";
  
  double massVal=startParams.first;
  double widthVal=startParams.second;
  
  double massErr=errParams.first;
  double widthErr=errParams.second;

  double massMin=massVal-4*massErr;
  if (massMin<0.) massMin=0.;

  double massMax=massVal+4*massErr;

  double widthMin=0.03;
  if (widthVal<0.04) widthMin=0.; 
  double widthMax=2*widthVal;
 
  upar.Add(massStr, massVal, massErr, massMin, massMax);
  upar.Add(widthStr, widthVal, widthErr, widthMin, widthMax);
}


void AbsPsi2STo2K2PiGamLh::setMnUsrParamsFlatteMass(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal, std::string key){

  if (key=="f980"){

  double start_f980M=startVal.Flatf980;
  double start_FlatgPiPi=startVal.Flatf980gPiPi;
  double start_FlatgKK=startVal.Flatf980gKK;
  double err_f980M=errVal.Flatf980;
  double err_FlatgPiPi=errVal.Flatf980gPiPi;
  double err_FlatgKK=errVal.Flatf980gKK;

  //now fill the fitParameterMap
  std::string massStr="f980Mass";
  std::string f980gPiPiStr="f980gPiPi";
  std::string f980gKKStr="f980gKK";

  double massMin=0.85;
  double massMax=1.2;
  upar.Add(massStr, start_f980M, err_f980M, massMin, massMax);

  upar.Add(f980gPiPiStr, start_FlatgPiPi, err_FlatgPiPi, 0., start_FlatgPiPi+3*err_FlatgPiPi);
  upar.Add(f980gKKStr, start_FlatgKK, err_FlatgKK, 0., start_FlatgKK+3*err_FlatgKK);
  }

  else { Alert << "Key: " << key << " not supported for setting up the MINUIT start mass parameters!!!" << endmsg;
    exit(1);
  }

}



void AbsPsi2STo2K2PiGamLh::checkFitParamVal(Psi2STo2K2PiGamData::fitParamVal& fitVal){
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> >  JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();

  for ( itJPCLS=JPCLSPsiToChiGam.begin(); itJPCLS!=JPCLSPsiToChiGam.end(); ++itJPCLS){
    std::pair<double, double> tmpParam=fitVal.PsiToChiGam[(*itJPCLS)];
    if ( tmpParam.first<-100. || tmpParam.second <-100.)
      {
	Alert << "Fit values for" << (*itJPCLS)->name() << "out of range" << endmsg;
      exit(1);
      }
  }
}

complex<double> AbsPsi2STo2K2PiGamLh::chiTo2K892Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K892, double K892Mass, double K892Width){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  for (Spin lamK892_1=-1; lamK892_1<=1; lamK892_1++){

    Spin lamK892_2=lamK892_1;
    
    complex<double> ampKpPi0= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi0[1][lamK892_1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    complex<double> ampKpPi1= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi1[1][lamK892_1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    complex<double> ampKmPi0= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi0[1][lamK892_2][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    complex<double> ampKmPi1= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi1[1][lamK892_2][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K892Mass, K892Width, 1));
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiTo2K892;
    for ( itChiTo2K892=ChiTo2K892.begin(); itChiTo2K892!=ChiTo2K892.end(); ++itChiTo2K892){

      boost::shared_ptr<const JPCLS> ChiTo2K892State=itChiTo2K892->first;      
      double theMag=itChiTo2K892->second.first;
      double thePhi=itChiTo2K892->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));

      complex<double> ChiTo2K892Tmpls=ampKpPi0*ampKmPi1;
      ChiTo2K892Tmpls+=ampKpPi1*ampKmPi0;      

      ChiTo2K892Tmpls=theMag*expiphi*ChiTo2K892Tmpls*sqrt(2*ChiTo2K892State->L+1)
	 *Clebsch(ChiTo2K892State->L,0,ChiTo2K892State->S, 0, 0, 0)
	 *Clebsch(1, lamK892_1, 1, -lamK892_2, ChiTo2K892State->S, 0);

      result+=ChiTo2K892Tmpls;		
    }
  }

  return result;
}






complex<double> AbsPsi2STo2K2PiGamLh::chiToK1K1Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK1K1, double K1Mass0, double K1Width0,  double K1Mass1, double K1Width1){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  for (Spin lamK1=-1; lamK1<=1; lamK1++){

    complex<double> amp_K10_KpPi0= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi0[1][lamK1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KpPi0= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi0[1][lamK1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));    


    complex<double> amp_K10_KpPi1= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi1[1][lamK1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KpPi1= 
      sqrt(3)*( conj(theData->DfKst1pToKpPi1[1][lamK1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));


    complex<double> amp_K10_KmPi0= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi0[1][lamK1][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KmPi0= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi0[1][lamK1][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));


    complex<double> amp_K10_KmPi1= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi1[1][lamK1][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1Mass0, K1Width0, 1));
    complex<double> amp_K11_KmPi1= 
      sqrt(3)*( conj(theData->DfKst1mToKmPi1[1][lamK1][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1Mass1, K1Width1, 1));
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
    for ( it=ChiToK1K1.begin(); it!=ChiToK1K1.end(); ++it){

      boost::shared_ptr<const JPCLS> ChiToK1K1State=it->first;      
      double theMag=it->second.first;
      double thePhi=it->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));

      complex<double> ChiToK1K1Tmpls=amp_K10_KpPi0*amp_K11_KmPi1;
      ChiToK1K1Tmpls+=amp_K10_KpPi1*amp_K11_KmPi0;      
      ChiToK1K1Tmpls+=amp_K10_KmPi0*amp_K11_KpPi1;
      ChiToK1K1Tmpls+=amp_K10_KmPi1*amp_K11_KpPi0;

      result+=theMag*expiphi*ChiToK1K1Tmpls*sqrt(2*ChiToK1K1State->L+1)
	 *Clebsch(ChiToK1K1State->L,0,ChiToK1K1State->S, 0, 0, 0)
	 *Clebsch(1, lamK1, 1, -lamK1, ChiToK1K1State->S, 0);

    }
  }

  return result;
}









complex<double> AbsPsi2STo2K2PiGamLh::chiTo2K_2_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_2, double K_2_Mass, double K_2_Width){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  for (Spin lamK2_1=-2; lamK2_1<=2; lamK2_1++){

    Spin lamK2_2=lamK2_1;
    
    complex<double> ampKpPi0= 
      sqrt(5)*( conj(theData->DfKst2pToKpPi0[2][lamK2_1][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    complex<double> ampKpPi1= 
      sqrt(5)*( conj(theData->DfKst2pToKpPi1[2][lamK2_1][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    complex<double> ampKmPi0= 
      sqrt(5)*( conj(theData->DfKst2mToKmPi0[2][lamK2_2][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    complex<double> ampKmPi1= 
      sqrt(5)*( conj(theData->DfKst2mToKmPi1[2][lamK2_2][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K_2_Mass, K_2_Width, 2));
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiTo2K2;
    for ( itChiTo2K2=ChiTo2K_2.begin(); itChiTo2K2!=ChiTo2K_2.end(); ++itChiTo2K2){

      boost::shared_ptr<const JPCLS> ChiTo2K_2_State=itChiTo2K2->first;      
      double theMag=itChiTo2K2->second.first;
      double thePhi=itChiTo2K2->second.second;
      complex<double> expiphi(cos(thePhi), sin(thePhi));

      complex<double> ChiTo2K2Tmpls=ampKpPi0*ampKmPi1;
      ChiTo2K2Tmpls+=ampKpPi1*ampKmPi0;      

      ChiTo2K2Tmpls=theMag*expiphi*ChiTo2K2Tmpls*sqrt(2*ChiTo2K_2_State->L+1)
	*Clebsch(ChiTo2K_2_State->L,0,ChiTo2K_2_State->S, 0, 0, 0)
	*Clebsch(2, lamK2_1, 2, -lamK2_2, ChiTo2K_2_State->S, 0);

      result+=ChiTo2K2Tmpls;		
    }
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK0K2Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK0K2, double K0Mass, double K0Width,  double K2Mass, double K2Width){
  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

    complex<double> ampK0ToKpPi0=  sqrt(3)*BreitWigner(KpPi0, K0Mass, K0Width);    
    complex<double> ampK0ToKpPi1=  sqrt(3)*BreitWigner(KpPi1, K0Mass, K0Width);
    complex<double> ampK0ToKmPi0=  sqrt(3)*BreitWigner(KmPi0, K0Mass, K0Width);    
    complex<double> ampK0ToKmPi1=  sqrt(3)*BreitWigner(KmPi1, K0Mass, K0Width);
    complex<double> ampK2ToKpPi0=  sqrt(5)*(conj(theData->DfKst2pToKpPi0[2][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K2Mass, K2Width, 2));
    complex<double> ampK2ToKpPi1=  sqrt(5)*(conj(theData->DfKst2pToKpPi1[2][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K2Mass, K2Width, 2));
    complex<double> ampK2ToKmPi0=  sqrt(5)*(conj(theData->DfKst2mToKmPi0[2][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K2Mass, K2Width, 2));
    complex<double> ampK2ToKmPi1=  sqrt(5)*(conj(theData->DfKst2mToKmPi1[2][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K2Mass, K2Width, 2));

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToK0K2.begin(); it!=ChiToK0K2.end(); ++it){
    
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));



    result+=theMag*expiphi*
      (ampK0ToKpPi0*ampK2ToKmPi1+ampK0ToKpPi1*ampK2ToKmPi0+ampK0ToKmPi0+ampK2ToKpPi1+ampK0ToKmPi1*ampK2ToKpPi0);
    
  } 
 
  return result;

}

complex<double> AbsPsi2STo2K2PiGamLh::chiTo2K_0_Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTo2K_0, double K_0_Mass0, double K_0_Width0, double K_0_Mass1, double K_0_Width1){

  complex<double> result(0.,0.);
  
  Vector4<double> KpPi0(theData->KpPi0_HeliChic0_4V.E(), theData->KpPi0_HeliChic0_4V.Px(), 
			 theData->KpPi0_HeliChic0_4V.Py(), theData->KpPi0_HeliChic0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliChic0_4V.E(), theData->KpPi1_HeliChic0_4V.Px(), 
			theData->KpPi1_HeliChic0_4V.Py(), theData->KpPi1_HeliChic0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliChic0_4V.E(), theData->KmPi0_HeliChic0_4V.Px(), 
			theData->KmPi0_HeliChic0_4V.Py(), theData->KmPi0_HeliChic0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliChic0_4V.E(), theData->KmPi1_HeliChic0_4V.Px(), 
			theData->KmPi1_HeliChic0_4V.Py(), theData->KmPi1_HeliChic0_4V.Pz());

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiTo2K0;
  for ( itChiTo2K0=ChiTo2K_0.begin(); itChiTo2K0!=ChiTo2K_0.end(); ++itChiTo2K0){
    
    double theMag=itChiTo2K0->second.first;
    double thePhi=itChiTo2K0->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));
    
    result+=theMag*expiphi*
      ( BreitWigner(KpPi0, K_0_Mass0, K_0_Width0) * BreitWigner(KmPi1, K_0_Mass1, K_0_Width1)
	+ BreitWigner(KpPi0, K_0_Mass1, K_0_Width1) * BreitWigner(KmPi1, K_0_Mass0, K_0_Width0)
	+ BreitWigner(KpPi1, K_0_Mass0, K_0_Width0) * BreitWigner(KmPi0, K_0_Mass1, K_0_Width1)
	+ BreitWigner(KpPi1, K_0_Mass1, K_0_Width1) * BreitWigner(KmPi0, K_0_Mass0, K_0_Width0) );
    
  } 
 
  return result;
}


complex<double> AbsPsi2STo2K2PiGamLh::chiToK1ToK1piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1ToK1dPi, double K1Mass, double K1Width, double K1dMass, double K1dWidth){
  complex<double> result(0.,0.);
  
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());

  Spin lamK1=0;
  Spin lamK1d=0;


  complex<double> resultTmpLoop(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itK1ToK1dPi;
  for ( itK1ToK1dPi=K1ToK1dPi.begin(); itK1ToK1dPi!=K1ToK1dPi.end(); ++itK1ToK1dPi){
    
    boost::shared_ptr<const JPCLS> K1ToK1PiState=itK1ToK1dPi->first;
    double theMag=itK1ToK1dPi->second.first;
    double thePhi=itK1ToK1dPi->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K1prodClebsch=sqrt(5.);         
    double K1decClebsch=sqrt(2.*K1ToK1PiState->L+1)*Clebsch(K1ToK1PiState->L,0, K1ToK1PiState->S, 0, K1ToK1PiState->J, 0)
      *sqrt(5.);  //K1d decay

    complex<double> ampKp1ToKp1dPi0=
      conj(theData->DfK1400pToKstpPi0[1][0][0])*BreitWignerBlattW(KpPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L)//K1dec
      *conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec

    complex<double> ampKp1ToKp1dPi1=
      conj(theData->DfK1400pToKstpPi1[1][0][0])*BreitWignerBlattW(KpPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L) //K1dec
      *conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
    
    complex<double> ampKm1ToKm1dPi0=
      conj(theData->DfK1400mToKstmPi0[1][0][0])*BreitWignerBlattW(KmPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L) //K1dec
      *conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec
    
    complex<double> ampKm1ToKm1dPi1=
      conj(theData->DfK1400mToKstmPi1[1][0][0])*BreitWignerBlattW(KmPiPi, K1dMass, 0.1349766, K1Mass, K1Width, K1ToK1PiState->L) //K1dec
      *conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

//     complex<double> ampKp1ToKp1dPi0=
//       conj(theData->DfK1400pToKstpPi0[1][0][0])*BreitWigner(KpPiPi, K1Mass, K1Width)//K1dec
//       *conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec

//     complex<double> ampKp1ToKp1dPi1=
//       conj(theData->DfK1400pToKstpPi1[1][0][0])*BreitWigner(KpPiPi, K1Mass, K1Width) //K1dec
//       *conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);
    
//     complex<double> ampKm1ToKm1dPi0=
//       conj(theData->DfK1400mToKstmPi0[1][0][0])*BreitWigner(KmPiPi, K1Mass, K1Width) //K1dec
//       *conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); //K1ddec
    
//     complex<double> ampKm1ToKm1dPi1=
//       conj(theData->DfK1400mToKstmPi1[1][0][0])*BreitWigner(KmPiPi, K1Mass, K1Width) //K1dec
//       *conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

    resultTmpLoop+=theMag*expiphi*K1prodClebsch*K1decClebsch*(ampKp1ToKp1dPi0+ampKp1ToKp1dPi1+ampKm1ToKm1dPi0+ampKm1ToKm1dPi1);
  }
  result=resultTmpLoop;
  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK1ToK0piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K1ToK0Pi, double K1Mass, double K1Width, double K0Mass, double K0Width){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itK1ToK0Pi;
  for ( itK1ToK0Pi=K1ToK0Pi.begin(); itK1ToK0Pi!=K1ToK0Pi.end(); ++itK1ToK0Pi){
    
    boost::shared_ptr<const JPCLS> K1ToK0PiState=itK1ToK0Pi->first;
    double theMag=itK1ToK0Pi->second.first;
    double thePhi=itK1ToK0Pi->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K1prodClebsch=sqrt(3.);         
    double K1decClebsch=sqrt(3.);  //K1d decay   

    complex<double> ampKpToKp0Pi0=
      conj(theData->DfK1400pToKstpPi0[1][0][0])*BreitWignerBlattW(KpPiPi, K0Mass, 0.1349766, K1Mass, K1Width, K1ToK0PiState->L)//K1dec
      *BreitWigner(KpPi1, K0Mass, K0Width); 

    complex<double> ampKp1400ToKp1dPi1=
      conj(theData->DfK1400pToKstpPi1[1][0][0])*BreitWignerBlattW(KpPiPi, K0Mass, 0.1349766, K1Mass, K1Width, K1ToK0PiState->L) //K1dec
      *BreitWigner(KpPi0, K0Mass, K0Width);
    
    complex<double> ampKm1400ToKm1dPi0=
      conj(theData->DfK1400mToKstmPi0[1][0][0])*BreitWignerBlattW(KmPiPi, K0Mass, 0.1349766, K1Mass, K1Width, K1ToK0PiState->L) //K1dec
      *BreitWigner(KmPi1, K0Mass, K0Width); //K1ddec
    
    complex<double> ampKm1400ToKp1dPi1=
      conj(theData->DfK1400mToKstmPi1[1][0][0])*BreitWignerBlattW(KmPiPi, K0Mass, 0.1349766, K1Mass, K1Width, K1ToK0PiState->L) //K1dec
      *BreitWigner(KmPi0, K0Mass, K0Width); //K1ddec

    result+=theMag*expiphi*K1prodClebsch*K1decClebsch*(ampKpToKp0Pi0+ampKp1400ToKp1dPi1+ampKm1400ToKm1dPi0+ampKm1400ToKp1dPi1);
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK0KT0KpipiKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K0dMass, double K0dWidth){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=K0ToK0Pi.begin(); it!=K0ToK0Pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> K0ToK0PiState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K1prodClebsch=sqrt(3.);         
    double K1decClebsch=sqrt(3.);  //K1d decay   

    complex<double> ampKpToKp0Pi0=BreitWigner(KpPiPi, K0pMass, K0pWidth)//K0prod
      *BreitWigner(KpPi1, K0dMass, K0dWidth); 

    complex<double> ampKpToKp1dPi1=BreitWigner(KpPiPi, K0pMass, K0pWidth) //K0prod
      *BreitWigner(KpPi0, K0dMass, K0dWidth);

    complex<double> ampKmToKm0Pi0=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *BreitWigner(KmPi1, K0dMass, K0dWidth);     

    complex<double> ampKmToKm1Pi1=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *BreitWigner(KmPi0, K0dMass, K0dWidth);

    result+=theMag*expiphi*K1prodClebsch*K1decClebsch*(ampKpToKp0Pi0+ampKpToKp1dPi1+ampKmToKm0Pi0+ampKmToKm1Pi1);
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK0KT0K1piKAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& K0ToK0Pi, double K0pMass, double K0pWidth, double K1dMass, double K1dWidth){
  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  Vector4<double> KpPi0(theData->KpPi0_HeliKpPi0Pi0_4V.E(), theData->KpPi0_HeliKpPi0Pi0_4V.Px(), 
			 theData->KpPi0_HeliKpPi0Pi0_4V.Py(), theData->KpPi0_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KpPi1(theData->KpPi1_HeliKpPi0Pi0_4V.E(), theData->KpPi1_HeliKpPi0Pi0_4V.Px(), 
			theData->KpPi1_HeliKpPi0Pi0_4V.Py(), theData->KpPi1_HeliKpPi0Pi0_4V.Pz());

  Vector4<double> KmPi0(theData->KmPi0_HeliKmPi0Pi0_4V.E(), theData->KmPi0_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi0_HeliKmPi0Pi0_4V.Py(), theData->KmPi0_HeliKmPi0Pi0_4V.Pz());

  Vector4<double> KmPi1(theData->KmPi1_HeliKmPi0Pi0_4V.E(), theData->KmPi1_HeliKmPi0Pi0_4V.Px(), 
			theData->KmPi1_HeliKmPi0Pi0_4V.Py(), theData->KmPi1_HeliKmPi0Pi0_4V.Pz());


  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=K0ToK0Pi.begin(); it!=K0ToK0Pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> K0ToK0PiState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    double K0proddecClebsch=2.*sqrt(3.);         
    double K1decClebsch=sqrt(5.);  //K1d decay   

    complex<double> ampKpToK1pPi0=BreitWigner(KpPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1pToKpPi1ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1); 


    complex<double> ampKpToK1pPi1=BreitWigner(KpPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1pToKpPi0ViaKpPiPi[1][0][0])*BreitWignerBlattW(KpPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);


    complex<double> ampKmToK1mPi0=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1mToKmPi1ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi1, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);

    complex<double> ampKmToK1mPi1=BreitWigner(KmPiPi, K0pMass, K0pWidth)//K0prod
      *conj(theData->DfKst1mToKmPi0ViaKmPiPi[1][0][0])*BreitWignerBlattW(KmPi0, 0.493677, 0.1349766, K1dMass, K1dWidth, 1);


    result+=theMag*expiphi*K0proddecClebsch*K1decClebsch*(ampKpToK1pPi0+ampKpToK1pPi1+ampKmToK1mPi0+ampKmToK1mPi1);
  }

  return result;
}



complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif0_k, double f980_piMass, double f980_piWidth, double f0_kMass, double f0_kWidth){

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
			 theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());

  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
			 theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_pif0_k.begin(); it!=ChiTof980_pif0_k.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(PiPi, f980_piMass,f980_piWidth)*BreitWigner(KK, f0_kMass, f0_kWidth);
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_pif0_k, double f980_piMass, double f980_g1, double f980_g2, double f0_kMass, double f0_kWidth){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_pif0_k.begin(); it!=ChiTof980_pif0_k.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_piMass, f980_g1, f980_g2)*BreitWigner(KK, f0_kMass, f0_kWidth);
  }

  return result;

}

complex<double> AbsPsi2STo2K2PiGamLh::chiTof980_kf0_piAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980_kf0_pi, double f980_Mass, double f980_gKK, double f980_gPiPi, double f0_piMass, double f0_piWidth){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;

  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980_kf0_pi.begin(); it!=ChiTof980_kf0_pi.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(KK, decPair2, decPair1, f980_Mass, f980_gKK, f980_gPiPi)*BreitWigner(PiPi, f0_piMass, f0_piWidth);
  }

  return result;

}



complex<double> AbsPsi2STo2K2PiGamLh::chiTof980f980Amp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof980980, double f980_piMass, double f980_g1, double f980_g2){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);

  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);

  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof980980.begin(); it!=ChiTof980980.end(); ++it){
    
    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_piMass, f980_g1, f980_g2)*Flatte(KK, decPair2, decPair1, f980_piMass, f980_g2, f980_g1);
  }

  return result;

}

complex<double>  AbsPsi2STo2K2PiGamLh::chiTof0_pif0_kAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiTof0_pif0_k, double f0_pi_Mass, double f0_pi_Width,  double f0_kMass, double f0_kWidth){
 
  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());
  
  Vector4<double> KK(theData->KpKm_HeliChic0_4V.E(), theData->KpKm_HeliChic0_4V.Px(), 
		     theData->KpKm_HeliChic0_4V.Py(), theData->KpKm_HeliChic0_4V.Pz());
  
  complex<double> result(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiTof0_pif0_k.begin(); it!=ChiTof0_pif0_k.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*BreitWigner(PiPi, f0_pi_Mass, f0_pi_Width)*BreitWigner(KK, f0_kMass, f0_kWidth);
  }

  return result;
}

complex<double> AbsPsi2STo2K2PiGamLh::chiToK_0_KToKf980KAmp(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& ChiToK_0_KToKf980K, double f980_Mass, double f980_gPiPi,  double f980_gKK, double K_0_Mass, double K_0_Width){

  const double massPi0 = 0.1349766;
  std::pair <const double, const double> decPair1=make_pair(massPi0, massPi0);
  
  const double KplusMass = 0.493677;
  const double K0Mass = 0.497614;
  
  std::pair <const double, const double> decPair2=make_pair(KplusMass,K0Mass);
  
  Vector4<double> PiPi(theData->PiPi_HeliChic0_4V.E(), theData->PiPi_HeliChic0_4V.Px(), 
		       theData->PiPi_HeliChic0_4V.Py(), theData->PiPi_HeliChic0_4V.Pz());


  Vector4<double> KpPiPi(theData->KpPiPi_HeliChic0_4V.E(), theData->KpPiPi_HeliChic0_4V.Px(), 
			 theData->KpPiPi_HeliChic0_4V.Py(), theData->KpPiPi_HeliChic0_4V.Pz());

  Vector4<double> KmPiPi(theData->KmPiPi_HeliChic0_4V.E(), theData->KmPiPi_HeliChic0_4V.Px(), 
			 theData->KmPiPi_HeliChic0_4V.Py(), theData->KmPiPi_HeliChic0_4V.Pz());

  complex<double> result(0.,0.);


  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator it;
  for ( it=ChiToK_0_KToKf980K.begin(); it!=ChiToK_0_KToKf980K.end(); ++it){

    boost::shared_ptr<const JPCLS> theState=it->first;
    double theMag=it->second.first;
    double thePhi=it->second.second;
    complex<double> expiphi(cos(thePhi), sin(thePhi));

    result+=theMag*expiphi*Flatte(PiPi, decPair1, decPair2, f980_Mass, f980_gPiPi, f980_gKK)*(BreitWigner(KpPiPi, K_0_Mass, K_0_Width)+BreitWigner(KmPiPi, K_0_Mass, K_0_Width));
  }

  return result;

}

void AbsPsi2STo2K2PiGamLh::print(std::ostream& os) const{
  os << "AbsPsi2STo2K2PiGamLh::print\n";
}

void AbsPsi2STo2K2PiGamLh::printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const{
  //  print fit paramss

  DebugMsg <<"\t PhaseSpace:" << theParamVal.phaseSpace << endmsg;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();

  for ( itJPCLS=JPCLSPsiToChiGam.begin(); itJPCLS!=JPCLSPsiToChiGam.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "ChiGam" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.PsiToChiGam[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }  

}

void AbsPsi2STo2K2PiGamLh::dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const{

  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }
std::string tmpString="PhaseSpace"+suffix;
  os << tmpString << "\t" << theParamVal.phaseSpace << "\t 0." << std::endl;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates(); 

 for ( itJPCLS=JPCLSPsiToChiGam.begin(); itJPCLS!=JPCLSPsiToChiGam.end(); ++itJPCLS){
   std::string tmpStringDec=(*itJPCLS)->name()+"ChiGam"+suffix;
    std::pair<double, double> tmpParam=theParamVal.PsiToChiGam[(*itJPCLS)];
    os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }
}
