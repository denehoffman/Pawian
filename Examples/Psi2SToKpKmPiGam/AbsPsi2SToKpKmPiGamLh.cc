#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsPsi2SToKpKmPiGamLh::AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList, boost::shared_ptr<const Psi2SToKpKmPiGamStates> theStates) :
  _Psi2SToKpKmPiGamEvtListPtr(theEvtList),
  _Psi2SToKpKmPiGamStatesPtr(theStates)
{
  _evtDataVec=_Psi2SToKpKmPiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2SToKpKmPiGamEvtListPtr->getMcVecs();
}

AbsPsi2SToKpKmPiGamLh::AbsPsi2SToKpKmPiGamLh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theAbsPsi2SToKpKmPiGamLhPtr):
  _Psi2SToKpKmPiGamEvtListPtr(theAbsPsi2SToKpKmPiGamLhPtr->getEventList()),
  _Psi2SToKpKmPiGamStatesPtr(theAbsPsi2SToKpKmPiGamLhPtr->getPsi2SToKpKmPiGamStates())
{
  _evtDataVec=_Psi2SToKpKmPiGamEvtListPtr->getDataVecs();
  _evtMCVec=_Psi2SToKpKmPiGamEvtListPtr->getMcVecs();
}

AbsPsi2SToKpKmPiGamLh::~AbsPsi2SToKpKmPiGamLh()
{
}

double AbsPsi2SToKpKmPiGamLh::calcLogLh(const Psi2SToKpKmPiGamData::fitParamVal& theParamVal){
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData*>::iterator iterm;
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

double AbsPsi2SToKpKmPiGamLh::calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal){

  Spin Psi2SM=1;
  Spin GamM=1;
  complex<double> AmpPsi2SMpGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=1;
  GamM=-1;
  complex<double> AmpPsi2SMpGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=-1;
  GamM=1; 
  complex<double> AmpPsi2SMmGp=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

  Psi2SM=-1;
  GamM=-1; 
  complex<double> AmpPsi2SMmGm=calcCoherentAmp(Psi2SM, GamM, theParamVal, theData);

//   DebugMsg << "AmpPsi2SMp " << AmpPsi2SMp << endmsg;
  
  double result=norm(AmpPsi2SMpGp)+norm(AmpPsi2SMpGm)+norm(AmpPsi2SMmGp)+norm(AmpPsi2SMmGm);

  return result;  
}



int AbsPsi2SToKpKmPiGamLh::setFitParamValDec(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const{

  int resultCount=counter;
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > currentStates;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >* currentMap=0;


  if (key=="ChiGam"){ 
    currentStates=_Psi2SToKpKmPiGamStatesPtr->PsiToChiGamStates();
    currentMap= &theParamVal.PsiToChiGam;
  }
  else if (key=="K890Pi"){
    currentStates=_Psi2SToKpKmPiGamStatesPtr->ChiToKst1KStates();
    currentMap= &theParamVal.ChiToK890K;
  }
  else if (key=="K_0_1400Pi"){
    currentStates=_Psi2SToKpKmPiGamStatesPtr->ChiToKst0KStates();
    currentMap= &theParamVal.ChiToK1400_0_K;
  }
  else if (key=="K_1_1400Pi"){
    currentStates=_Psi2SToKpKmPiGamStatesPtr->ChiToKst1KStates();
    currentMap= &theParamVal.ChiToK1400_1_K;
  }
  else if (key=="K_2_1400Pi"){
    currentStates=_Psi2SToKpKmPiGamStatesPtr->ChiToKst2KStates();
    currentMap= &theParamVal.ChiToK1400_2_K;
  }
  else if (key=="a980Pi"){
    currentStates=_Psi2SToKpKmPiGamStatesPtr->ChiToa0PiStates();
    currentMap= &theParamVal.ChiToa0Pi;
  }
  else if (key=="a2Pi"){
    currentStates=_Psi2SToKpKmPiGamStatesPtr->ChiToa2PiStates();
    currentMap= &theParamVal.ChiToa2Pi;
  }
  else { Alert << "Key: " << key << " not supported for setting up the fit amplitude parameters!!!" << endmsg;
    exit(1);
  }

  for ( itJPCLS=currentStates.begin(); itJPCLS!=currentStates.end(); ++itJPCLS){
    double mag=par[resultCount];
    resultCount++;
    double phi=0.;
    if (resultCount>1 || key != "ChiGam"){ phi=par[resultCount];
      resultCount++;
    }
    std::pair <double,double> tmpParameter=make_pair(mag,phi);
    (*currentMap)[(*itJPCLS)]=tmpParameter;
  }

  return resultCount;
}


int AbsPsi2SToKpKmPiGamLh::setFitParamValMass(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const{

  int resultCount=counter;

  pair<double, double>* currentPair=0;

  if (key=="K890"){
    currentPair= &theParamVal.BwK890;
  }
  else if (key=="K_0_1400"){
   currentPair= &theParamVal.BwK1400_0;
  }
  else if (key=="K_1_1400"){
    currentPair= &theParamVal.BwK1400_1;
  }
  else if (key=="K_2_1400"){
    currentPair= &theParamVal.BwK1400_2;
  }
  else if (key=="a980"){
    currentPair= &theParamVal.Bwa980;
  }
  else if (key=="a2"){
    currentPair= &theParamVal.Bwa2;
  }
  else { Alert << "Key: " << key << " not supported for setting up the fit mass parameters!!!" << endmsg;
    exit(1);
  }

  double currentMass=par[resultCount];
  resultCount++;
  double currentWidth=par[resultCount];
  resultCount++;

  std::pair <double,double> thePair=make_pair(currentMass, currentWidth);

  (*currentPair)=thePair;

  return resultCount;
}

int AbsPsi2SToKpKmPiGamLh::setFitParamFlattea980Mass(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par, int counter, std::string key) const{

  double* resultFlatMa980=0;
  double* resultFlatgKK=0;
  double* resultFlatgEtaPi=0;
  int resultCount=counter;

  if (key=="Flattea980"){
    resultFlatMa980= &theParamVal.FlatMa980;
    resultFlatgKK= &theParamVal.FlatgKK;
    resultFlatgEtaPi= &theParamVal.FlatgEtaPi;
  }
  else { Alert << "Key: " << key << " not supported for setting up the Flatte mass parameters!!!" << endmsg;
    exit(1);
  }

  double currentMass=par[resultCount];
  resultCount++;
  (*resultFlatMa980)=currentMass;

  double currentFlatgKK=par[resultCount];
  resultCount++;
  (*resultFlatgKK)=currentFlatgKK;

  double currentFlatgEtaPi=par[resultCount];
  resultCount++;
  (*resultFlatgEtaPi)=currentFlatgEtaPi;

  return resultCount;
}



complex<double> AbsPsi2SToKpKmPiGamLh::a980Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi, double a980Mass, double a980Width, Spin& lamChi){

  Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());
  

  complex<double> result=conj(theData->DfChiToa0Pi[1][lamChi][0])*BreitWigner(KpKm4V, a980Mass, a980Width); // Amp(a0(980)->K+K-)=1;
  
  complex<double> ChiToA980pi0AmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToa0Pi;
  
  for ( itChiToa0Pi=ChiToa0Pi.begin(); itChiToa0Pi!=ChiToa0Pi.end(); ++itChiToa0Pi){
    boost::shared_ptr<const JPCLS> ChiToa0State=itChiToa0Pi->first;
    double theChiToa0Mag=itChiToa0Pi->second.first;
    double theChiToa0Phi=itChiToa0Pi->second.second;
    complex<double> expiphiChiToa0(cos(theChiToa0Phi), sin(theChiToa0Phi));
    
    ChiToA980pi0AmpTmpls+=theChiToa0Mag*expiphiChiToa0*sqrt(2*ChiToa0State->L+1); 
    //Clebsch(ChiToa0State->L=1,0,ChiToa0State->S=0, a0lam-pilam=0, ChiToa0State->J=1, 0)=1  
    //Clebsch(a0.J=0, a0lam=0, piJ=0, pilam=0, ChiToa0State->S=0, a0lam-pilam=0)=1
  }
  
  result*=ChiToA980pi0AmpTmpls;

  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::a980FlatteAmp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi, double a980FlatteMass, double a980FlatgKK, double a980FlatgEtaPi, Spin& lamChi){

  Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());
  

//   complex<double> result=conj(theData->DfChiToa0Pi[1][lamChi][0])*BreitWigner(KpKm4V, a980Mass, a980Width); // Amp(a0(980)->K+K-)=1;

  complex<double> result=conj(theData->DfChiToa0Pi[1][lamChi][0])*FlatteA980(KpKm4V, a980FlatteMass, a980FlatgKK, a980FlatgEtaPi); // Amp(a0(980)->K+K-)=1;  

  complex<double> ChiToA980pi0AmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToa0Pi;
  
  for ( itChiToa0Pi=ChiToa0Pi.begin(); itChiToa0Pi!=ChiToa0Pi.end(); ++itChiToa0Pi){
    boost::shared_ptr<const JPCLS> ChiToa0State=itChiToa0Pi->first;
    double theChiToa0Mag=itChiToa0Pi->second.first;
    double theChiToa0Phi=itChiToa0Pi->second.second;
    complex<double> expiphiChiToa0(cos(theChiToa0Phi), sin(theChiToa0Phi));
    
    ChiToA980pi0AmpTmpls+=theChiToa0Mag*expiphiChiToa0*sqrt(2*ChiToa0State->L+1); 
    //Clebsch(ChiToa0State->L=1,0,ChiToa0State->S=0, a0lam-pilam=0, ChiToa0State->J=1, 0)=1  
    //Clebsch(a0.J=0, a0lam=0, piJ=0, pilam=0, ChiToa0State->S=0, a0lam-pilam=0)=1
  }
  
  result*=ChiToA980pi0AmpTmpls;

  return result;
}


complex<double> AbsPsi2SToKpKmPiGamLh::a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());

  for (Spin lam_a2=-2; lam_a2<=2; lam_a2++){
    
    complex<double> ChiToa2PiAmpTmp1=conj(theData->DfChiToa2Pi[1][lamChi][lam_a2])
      *sqrt(5)*conj(theData->Dfa2[2][lam_a2][0])*BreitWigner(KpKm4V, a2Mass, a2Width);
    
    complex<double> ChiToa2PiAmpTmpls(0.,0.);
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToa2Pi;
    for ( itChiToa2Pi=ChiToa2Pi.begin(); itChiToa2Pi!=ChiToa2Pi.end(); ++itChiToa2Pi){
	    boost::shared_ptr<const JPCLS> ChiToa2PiState=itChiToa2Pi->first;
	    
	    if ( fabs(lam_a2)>ChiToa2PiState->J || fabs(lam_a2)>ChiToa2PiState->S) continue;
	    
	    double theChiToa2PiMag=itChiToa2Pi->second.first;
	    double theChiToa2PiPhi=itChiToa2Pi->second.second;
	    complex<double> expiphiChiToa2Pi(cos(theChiToa2PiPhi), sin(theChiToa2PiPhi));
	    
	    
	    ChiToa2PiAmpTmpls+=theChiToa2PiMag*expiphiChiToa2Pi*sqrt(2*ChiToa2PiState->L+1)
	      *Clebsch(ChiToa2PiState->L,0,ChiToa2PiState->S, lam_a2, ChiToa2PiState->J, lam_a2)
	      *Clebsch(2, lam_a2, 0, 0, ChiToa2PiState->S, lam_a2);
          }
	  
	  ChiToa2PiAmpTmp1*=ChiToa2PiAmpTmpls;
	  result+=ChiToa2PiAmpTmp1;
  }

  return result;
}



complex<double> AbsPsi2SToKpKmPiGamLh::K892Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K, double K892Mass, double K892Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

  for (Spin lamK890=-1; lamK890<=1; lamK890++){
    
    complex<double> ChiToK890KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK890])
      *sqrt(3)*conj(theData->DfKp1[1][lamK890][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K892Mass, K892Width, 1)
      +conj(theData->DfChiToKmPi[1][lamChi][lamK890])
      *sqrt(3)*conj(theData->DfKm1[1][lamK890][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K892Mass, K892Width, 1);	
    
    complex<double> ChiToK890KAmpTmpls(0.,0.);
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToK890K;
    
    for ( itChiToK890K=ChiToK890K.begin(); itChiToK890K!=ChiToK890K.end(); ++itChiToK890K){
      boost::shared_ptr<const JPCLS> ChiToK890KState=itChiToK890K->first;
      
      if (lamK890>ChiToK890KState->J || lamK890 > ChiToK890KState->S) continue;
      
      double theChiToK890KMag=itChiToK890K->second.first;
      double theChiToK890KPhi=itChiToK890K->second.second;
      complex<double> expiphiChiToK890K(cos(theChiToK890KPhi), sin(theChiToK890KPhi));
      
      ChiToK890KAmpTmpls+=theChiToK890KMag*expiphiChiToK890K*sqrt(2*ChiToK890KState->L+1)
	*Clebsch(ChiToK890KState->L,0,ChiToK890KState->S, lamK890, ChiToK890KState->J, lamK890)
	*Clebsch(1, lamK890, 0, 0, ChiToK890KState->S, lamK890);
    }
    
    ChiToK890KAmpTmp1*=ChiToK890KAmpTmpls;
    result+=ChiToK890KAmpTmp1;
  }
  
  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::K0_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_0_K, double K0_1400Mass, double K0_1400Width, Spin& lamChi){

  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

  complex<double> result=conj(theData->DfChiToKpPi[1][lamChi][0])
    *BreitWigner(KpPi4V, K0_1400Mass, K0_1400Width)
    +conj(theData->DfChiToKmPi[1][lamChi][0])
    *BreitWigner(KmPi4V, K0_1400Mass, K0_1400Width);

  complex<double> ChiToK1400_0_KAmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToK1400_0_K;
  for ( itChiToK1400_0_K=ChiToK1400_0_K.begin(); itChiToK1400_0_K!=ChiToK1400_0_K.end(); ++itChiToK1400_0_K){
    boost::shared_ptr<const JPCLS> ChiToK1400_0_KState=itChiToK1400_0_K->first;
    double theChiToK1400_0_KMag=itChiToK1400_0_K->second.first;
    double theChiToK1400_0_KPhi=itChiToK1400_0_K->second.second;
    complex<double> expiphiChiToK1400_0_K(cos(theChiToK1400_0_KPhi), sin(theChiToK1400_0_KPhi));
    
    ChiToK1400_0_KAmpTmpls+=theChiToK1400_0_KMag*expiphiChiToK1400_0_K*sqrt(2*ChiToK1400_0_KState->L+1);
    //      *Clebsch(ChiToK1400_0_KState->L,0,ChiToK1400_0_KState->S, 0, ChiToK1400_0_KState->J, 0)
    //      *Clebsch(0, 0, 0, 0, 0, 0)
  }

  result*=ChiToK1400_0_KAmpTmpls;

  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::K1_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1_1400_K, double K1_1400Mass, double K1_1400Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

  for (Spin lamK1400=-1; lamK1400<=1; lamK1400++){
    
//     complex<double> ChiToK1_1400KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1400])
//       *sqrt(3)*conj(theData->DfKp1[1][lamK1400][0])*BreitWigner(KpPi4V, K1_1400Mass, K1_1400Width)
//       +conj(theData->DfChiToKmPi[1][lamChi][lamK1400])
//       *sqrt(3)*conj(theData->DfKm1[1][lamK1400][0])*BreitWigner(KmPi4V, K1_1400Mass, K1_1400Width);

    complex<double> ChiToK1_1400KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1400])
      *sqrt(3)*conj(theData->DfKp1[1][lamK1400][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K1_1400Mass, K1_1400Width, 1)
      +conj(theData->DfChiToKmPi[1][lamChi][lamK1400])
      *sqrt(3)*conj(theData->DfKm1[1][lamK1400][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K1_1400Mass, K1_1400Width, 1);
    complex<double> ChiToK1_1400KAmpTmpls(0.,0.);

    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToK1_1400K;
    for ( itChiToK1_1400K=ChiToK1_1400_K.begin(); itChiToK1_1400K!=ChiToK1_1400_K.end(); ++itChiToK1_1400K){
      boost::shared_ptr<const JPCLS> ChiToK1_1400KState=itChiToK1_1400K->first;
      double theChiToK1_1400KMag=itChiToK1_1400K->second.first;
      double theChiToK1_1400KPhi=itChiToK1_1400K->second.second;
      complex<double> expiphiChiToK1_1400K(cos(theChiToK1_1400KPhi), sin(theChiToK1_1400KPhi));
      
      ChiToK1_1400KAmpTmpls+=theChiToK1_1400KMag*expiphiChiToK1_1400K*sqrt(2*ChiToK1_1400KState->L+1)
	*Clebsch(ChiToK1_1400KState->L,0,ChiToK1_1400KState->S, lamK1400, ChiToK1_1400KState->J, lamK1400)
	*Clebsch(1, lamK1400, 0, 0, ChiToK1_1400KState->S, lamK1400);
    }
    
    ChiToK1_1400KAmpTmp1*=ChiToK1_1400KAmpTmpls;
    result+=ChiToK1_1400KAmpTmp1;
  }

  return result;
}

complex<double> AbsPsi2SToKpKmPiGamLh::K2_1400Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_2_K, double K2_1400Mass, double K2_1400Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V(theData->KpPi_HeliChic1_4V.E(), theData->KpPi_HeliChic1_4V.Px(), theData->KpPi_HeliChic1_4V.Py(),  theData->KpPi_HeliChic1_4V.Pz());
  Vector4<double> KmPi4V(theData->KmPi_HeliChic1_4V.E(), theData->KmPi_HeliChic1_4V.Px(), theData->KmPi_HeliChic1_4V.Py(),  theData->KmPi_HeliChic1_4V.Pz());

  for (Spin lamK1400_2=-2; lamK1400_2<=2; lamK1400_2++){
    
    complex<double> ChiToK1400_2_KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1400_2])
      *sqrt(5)*conj(theData->DfKp2[2][lamK1400_2][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K2_1400Mass, K2_1400Width, 2)
      +conj(theData->DfChiToKmPi[1][lamChi][lamK1400_2])
      *sqrt(5)*conj(theData->DfKm2[2][lamK1400_2][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K2_1400Mass, K2_1400Width, 2);
    
    complex<double> ChiToK1400_2_KAmpTmpls(0.,0.);
    
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToK1400_2_K;
    for ( itChiToK1400_2_K=ChiToK1400_2_K.begin(); itChiToK1400_2_K!=ChiToK1400_2_K.end(); ++itChiToK1400_2_K){
	    boost::shared_ptr<const JPCLS> ChiToK1400_2_KState=itChiToK1400_2_K->first;
	    
	    if ( fabs(lamK1400_2)>ChiToK1400_2_KState->J || fabs(lamK1400_2)>ChiToK1400_2_KState->S) continue;
	    
	    double theChiToK1400_2_KMag=itChiToK1400_2_K->second.first;
	    double theChiToK1400_2_KPhi=itChiToK1400_2_K->second.second;
	    complex<double> expiphiChiToK1400_2_K(cos(theChiToK1400_2_KPhi), sin(theChiToK1400_2_KPhi));
	    
	    
	    ChiToK1400_2_KAmpTmpls+=theChiToK1400_2_KMag*expiphiChiToK1400_2_K*sqrt(2*ChiToK1400_2_KState->L+1)
	      *Clebsch(ChiToK1400_2_KState->L,0,ChiToK1400_2_KState->S, lamK1400_2, ChiToK1400_2_KState->J, lamK1400_2)
	      *Clebsch(2, lamK1400_2, 0, 0, ChiToK1400_2_KState->S, lamK1400_2);
          }
	  
	  ChiToK1400_2_KAmpTmp1*=ChiToK1400_2_KAmpTmpls;
	  result+=ChiToK1400_2_KAmpTmp1;
  }

  return result;
}

void AbsPsi2SToKpKmPiGamLh::setMnUsrParamsDec(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal, std::string key){
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > startParams;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > errParams;

  if (key=="ChiGam"){
    startParams=startVal.PsiToChiGam;
    errParams=errVal.PsiToChiGam;
  }
  else if (key=="K890Pi"){
    startParams=startVal.ChiToK890K;
    errParams=errVal.ChiToK890K;
  }
  else if (key=="K_0_1400Pi"){
    startParams=startVal.ChiToK1400_0_K;
    errParams=errVal.ChiToK1400_0_K;
  }
  else if (key=="K_1_1400Pi"){
    startParams=startVal.ChiToK1400_1_K;
    errParams=errVal.ChiToK1400_1_K;
  }
  else if (key=="K_2_1400Pi"){
    startParams=startVal.ChiToK1400_2_K;
    errParams=errVal.ChiToK1400_2_K;
  }
  else if (key=="a980Pi"){
    startParams=startVal.ChiToa0Pi;
    errParams=errVal.ChiToa0Pi;
  }
  else if (key=="a2Pi"){
    startParams=startVal.ChiToa2Pi;
    errParams=errVal.ChiToa2Pi;
  }
  else { Alert << "Key: " << key << " not supported for setting up the MINUIT start amplitude parameters!!!" << endmsg;
    exit(1);
  }

  int counter=0;
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
    if (counter>0 || key != "ChiGam") upar.Add(phiStr, phiVal, phiErr, -2*M_PI, 2*M_PI);

    counter++;
  }

}

void AbsPsi2SToKpKmPiGamLh::setMnUsrParamsMass(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal, std::string key){
  pair<double, double> startParams;
  pair<double, double> errParams;

  if (key=="K890"){
    startParams=startVal.BwK890;
    errParams=errVal.BwK890;
  }
  else if (key=="K_0_1400"){
    startParams=startVal.BwK1400_0;
    errParams=errVal.BwK1400_0;
  }
  else if (key=="K_1_1400"){
    startParams=startVal.BwK1400_1;
    errParams=errVal.BwK1400_1;
  }
  else if (key=="K_2_1400"){
    startParams=startVal.BwK1400_2;
    errParams=errVal.BwK1400_2;
  }
  else if (key=="a980"){
    startParams=startVal.Bwa980;
    errParams=errVal.Bwa980;
  }
  else if (key=="a2"){
    startParams=startVal.Bwa2;
    errParams=errVal.Bwa2;
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

  double massMin=massVal-widthVal;
  if (massMin<0.) massMin=0.;

  double massMax=massVal+widthVal;

  double widthMin=0.; 
  double widthMax=2*widthVal;
 
  upar.Add(massStr, massVal, massErr, massMin, massMax);
  upar.Add(widthStr, widthVal, widthErr, widthMin, widthMax);
}


void AbsPsi2SToKpKmPiGamLh::setMnUsrParamsFlattea980Mass(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal, std::string key){
  double start_a980M;
  double start_FlatgKK;
  double start_FlatgEtaPi;
  double err_a980M;
  double err_FlatgKK;
  double err_FlatgEtaPi;

  start_a980M=startVal.FlatMa980;
  start_FlatgKK=startVal.FlatgKK;
  start_FlatgEtaPi=startVal.FlatgEtaPi;

  err_a980M=errVal.FlatMa980;
  err_FlatgKK=errVal.FlatgKK;
  err_FlatgEtaPi=errVal.FlatgEtaPi;
 
  //now fill the fitParameterMap
  std::string massStr="a980Mass";
  std::string gKKStr="gKK";
  std::string gEtaPiStr="gEtaPi";

  double massMin=0.97;
  double massMax=1.2;
  upar.Add(massStr, start_a980M, err_a980M, massMin, massMax);

  upar.Add(gKKStr, start_FlatgKK, err_FlatgKK, 0., start_FlatgKK+3*err_FlatgKK);
  upar.Add(gEtaPiStr, start_FlatgEtaPi, err_FlatgEtaPi, 0., start_FlatgEtaPi+3*err_FlatgEtaPi);
}



void AbsPsi2SToKpKmPiGamLh::checkFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& fitVal){
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> >  JPCLSPsiToChiGam=_Psi2SToKpKmPiGamStatesPtr->PsiToChiGamStates();

  for ( itJPCLS=JPCLSPsiToChiGam.begin(); itJPCLS!=JPCLSPsiToChiGam.end(); ++itJPCLS){
    std::pair<double, double> tmpParam=fitVal.PsiToChiGam[(*itJPCLS)];
    if ( tmpParam.first<-100. || tmpParam.second <-100.)
      {
	Alert << "Fit values for" << (*itJPCLS)->name() << "out of range" << endmsg;
      exit(1);
      }
  }
}

void AbsPsi2SToKpKmPiGamLh::print(std::ostream& os) const{
  os << "AbsPsi2SToKpKmPiGamLh::print\n";
}

void AbsPsi2SToKpKmPiGamLh::printCurrentFitResult(Psi2SToKpKmPiGamData::fitParamVal& theParamVal) const{
  //  print fit paramss
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2SToKpKmPiGamStatesPtr->PsiToChiGamStates();

  for ( itJPCLS=JPCLSPsiToChiGam.begin(); itJPCLS!=JPCLSPsiToChiGam.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "ChiGam" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.PsiToChiGam[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }  

}
