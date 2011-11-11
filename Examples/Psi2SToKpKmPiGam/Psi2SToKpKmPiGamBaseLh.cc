#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamBaseLh.hh"
//#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamBaseLh::Psi2SToKpKmPiGamBaseLh(boost::shared_ptr<const EvtDataBaseList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsLh(theEvtList)
  ,_K0_1430Hyp(true)
  ,_K1_1410Hyp(true)
  ,_K2_1430Hyp(true)
  ,_K1_1680Hyp(false)
  ,_KKPi_Hyp(false)
{

  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K0_1430HypBase");

  if (iter !=hypMap.end()){
    _K0_1430Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K0_1430Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K0_1430HypBase not set!!!" <<endmsg;

  iter= hypMap.find("K1_1410HypBase");
  if (iter !=hypMap.end()){
    _K1_1410Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1410Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1410HypBase not set!!!" <<endmsg;

  iter= hypMap.find("K2_1430HypBase");
  if (iter !=hypMap.end()){
    _K2_1430Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K2_1430Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K2_1430HypBase not set!!!" <<endmsg;

  iter= hypMap.find("KKPi_HypBase");
  if (iter !=hypMap.end()){
    _KKPi_Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _KKPi_Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis KKPi_HypBase not set!!!" <<endmsg;

  iter= hypMap.find("K1_1680_HypBase");
  if (iter !=hypMap.end()){
    _K1_1680Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1680Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }

  else Alert << "hypothesis KKPi_HypBase not set!!!" <<endmsg;
}

Psi2SToKpKmPiGamBaseLh::Psi2SToKpKmPiGamBaseLh( boost::shared_ptr<AbsLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsLh(theLhPtr->getEventList())
  ,_K0_1430Hyp(true)
  ,_K1_1410Hyp(true)
  ,_K2_1430Hyp(true)
  ,_K1_1680Hyp(false)
  ,_KKPi_Hyp(false)
{
  std::map<const std::string, bool>::const_iterator iter= hypMap.find("K0_1430HypBase");

  if (iter !=hypMap.end()){
    _K0_1430Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K0_1430Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K0_1430HypBase not set!!!" <<endmsg;

  iter= hypMap.find("K1_1410HypBase");
  if (iter !=hypMap.end()){
    _K1_1410Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1410Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K1_1410HypBase not set!!!" <<endmsg;

  iter= hypMap.find("K2_1430HypBase");
  if (iter !=hypMap.end()){
    _K2_1430Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K2_1430Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis K2_1430HypBase not set!!!" <<endmsg;

  iter= hypMap.find("KKPi_HypBase");
  if (iter !=hypMap.end()){
    _KKPi_Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _KKPi_Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else Alert << "hypothesis KKPi_HypBase not set!!!" <<endmsg;

  iter= hypMap.find("K1_1680_HypBase");
  if (iter !=hypMap.end()){
    _K1_1680Hyp= iter->second;
    Info<< "hypothesis " << iter->first << "\t" << _K1_1680Hyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
}

Psi2SToKpKmPiGamBaseLh::~Psi2SToKpKmPiGamBaseLh()
{;
}

double Psi2SToKpKmPiGamBaseLh::calcEvtIntensity(EvtData* theData, fitParams& theParamVal){


  double phaseSpaceVal=theParamVal.otherParams[paramEnumChic1ToKpKmPiGam::phaseSpace];  


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

//    DebugMsg << "AmpPsi2SMpGp " << AmpPsi2SMpGp << endmsg;
  
  double result=norm(AmpPsi2SMpGp)+norm(AmpPsi2SMpGm)+norm(AmpPsi2SMmGp)+norm(AmpPsi2SMmGm)+phaseSpaceVal;
  return result;  
}



complex<double> Psi2SToKpKmPiGamBaseLh::calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, EvtData* theData){

  complex<double> result(0.,0.);

      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& thePsiMags=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::ChiGam];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& thePsiPhis=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::ChiGam]; 

  for (Spin lamChi=-1; lamChi<=1; lamChi++){

      Spin lamChiGam=lamChi-lamGam;

      complex<double> PsiDecAmp=conj(theData->WignerDs[enumChic1ToKpKmPiGamData::dfPsi][1][Minit][lamChiGam]);
      complex<double> PsiDecAmpTmpls(0.,0.);  
      
 
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itPsi;      
      for ( itPsi=thePsiMags.begin(); itPsi!=thePsiMags.end(); ++itPsi){
    
	boost::shared_ptr<const JPCLS> PsiState=itPsi->first;

        if (fabs(lamChiGam)>PsiState->J || fabs(lamChiGam)>PsiState->S) continue;
	
	double thePsiMag=itPsi->second;
	double thePsiPhi=thePsiPhis[PsiState];
	complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
    
	complex<double> dummyAmp(0.,0.);
	complex<double> decAmp=calcDecAmp(dummyAmp, lamChi, theParamVal, theData);
	PsiDecAmpTmpls+=thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
	  *Clebsch(PsiState->L,0,PsiState->S,lamChiGam,PsiState->J,lamChiGam)
 	  *Clebsch(1, lamChi, 1, -lamGam, PsiState->S,lamChiGam)
	  *decAmp;
      }

      PsiDecAmp*=PsiDecAmpTmpls;
      result+=PsiDecAmp;
  }
      return result;
}

complex<double> Psi2SToKpKmPiGamBaseLh::calcDecAmp(complex<double>& inAmp,Spin lamChi, fitParams& theParamVal, EvtData* theData){

  complex<double> result=inAmp;

  //Chi to a0(980) pi0 and a0(980)->K+K-
  result+=a980FlatteAmp(theData, theParamVal, lamChi);
  
  //Chi -> K890 K and K890->K pi0
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK890KMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::K890K];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK890KPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::K890K];
  double K892Mass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::K890];
  double K892Width=theParamVal.Widths[paramEnumChic1ToKpKmPiGam::K890];
  
  //  result+=K892AmpK0_1400Amp(theData, ChiToK890KMag, ChiToK890KPhi, K892Mass, K892Width, lamChi);
  result+=K1_Amp(theData, ChiToK890KMag, ChiToK890KPhi, K892Mass, K892Width, lamChi);  
  //Chi -> K_0_1400 K and K_0_1400->K pi0
  if(_K0_1430Hyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_1400KMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::K_0_1400K];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_1400KPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::K_0_1400K];
    double K_0_1400Mass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::K_0_1400];
    double K_0_1400Width=theParamVal.Widths[paramEnumChic1ToKpKmPiGam::K_0_1400];
    
    result+=K0_Amp(theData, ChiToK_0_1400KMag, ChiToK_0_1400KPhi, K_0_1400Mass, K_0_1400Width, lamChi);
  }
  
  
  //Chi -> K_1_1400 K and K_1_1400->K pi0
  if(_K1_1410Hyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_1_1400KMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::K_1_1400K];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK__1400KPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::K_1_1400K];
    double K_1_1400Mass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::K_1_1400];
    double K_1_1400Width=theParamVal.Widths[paramEnumChic1ToKpKmPiGam::K_1_1400];
    
    result+=K1_Amp(theData, ChiToK_1_1400KMag, ChiToK__1400KPhi, K_1_1400Mass, K_1_1400Width, lamChi);
  }
  
  //Chi -> K_2_1400 K and K_2_1400->K pi0
  if(_K2_1430Hyp){
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_2_1400KMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::K_2_1400K];
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_2_1400KPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::K_2_1400K];
    double K_2_1400Mass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::K_2_1400];
    double K_2_1400Width=theParamVal.Widths[paramEnumChic1ToKpKmPiGam::K_2_1400];
    result+=K2_Amp(theData, ChiToK_2_1400KMag,  ChiToK_2_1400KPhi, K_2_1400Mass, K_2_1400Width, lamChi);
  }
  

//   //K_1_1680->K pi0
//   if(_K1_1680Hyp){
//     std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1680K=theParamVal.ChiToK1680_1_K; 
//     double K_1_1680Mass=theParamVal.BwK1680_1.first;
//     double K_1_1680Width=theParamVal.BwK1680_1.second;
//     result+=K892Amp(theData, ChiToK_1_1680K, K_1_1680Mass, K_1_1680Width, lamChi);
//   }

//   //Chi -> K+ K- pi0
//   if(_KKPi_Hyp){
//     std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi=theParamVal.ChiToKKPi;
//     result+=KKPi_Amp(theData, ChiToKKPi, lamChi);
//   }
  
  return result;
}


void Psi2SToKpKmPiGamBaseLh::print(std::ostream& os) const{
  os << "Psi2SToKpKmPiGamBaseLh::print\n";
}

complex<double> Psi2SToKpKmPiGamBaseLh::a980Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiPhi, double a980Mass, double a980Width, Spin& lamChi){

  Vector4<double>& KpKm4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1];

  complex<double> result=conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToa0Pi][1][lamChi][0])*BreitWigner(KpKm4V, a980Mass, a980Width); // Amp(a0(980)->K+K-)=1;
  
  complex<double> ChiToA980pi0AmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToa0PiMag;
  
  for ( itChiToa0PiMag=ChiToa0PiMag.begin(); itChiToa0PiMag!=ChiToa0PiMag.end(); ++itChiToa0PiMag){
    boost::shared_ptr<const JPCLS> ChiToa0State=itChiToa0PiMag->first;
    double theChiToa0Mag=itChiToa0PiMag->second;
    double theChiToa0Phi=ChiToa0PiPhi[ChiToa0State];
    complex<double> expiphiChiToa0(cos(theChiToa0Phi), sin(theChiToa0Phi));
    
    ChiToA980pi0AmpTmpls+=theChiToa0Mag*expiphiChiToa0*sqrt(2*ChiToa0State->L+1); 
    //Clebsch(ChiToa0State->L=1,0,ChiToa0State->S=0, a0lam-pilam=0, ChiToa0State->J=1, 0)=1
    //Clebsch(a0.J=0, a0lam=0, piJ=0, pilam=0, ChiToa0State->S=0, a0lam-pilam=0)=1
  }
  
  result*=ChiToA980pi0AmpTmpls;

  return result;
}


complex<double> Psi2SToKpKmPiGamBaseLh::a980FlatteAmp(EvtData* theData, fitParams& theParamVal, Spin& lamChi){

  Vector4<double>& KpKm4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpKm_HeliChic1];

  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToa0PiMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::a980Pi];
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToa0PiPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::a980Pi];  
  double a980FlatteMass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::a980];
  double a980FlatgKK=theParamVal.gFactors[paramEnumChic1ToKpKmPiGam::a980gKK];
  double a980FlatgEtaPi=theParamVal.gFactors[paramEnumChic1ToKpKmPiGam::a980gEtaPi];

  complex<double> result=conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToa0Pi][1][lamChi][0])*FlatteA980(KpKm4V, a980FlatteMass, a980FlatgKK, a980FlatgEtaPi); // Amp(a0(980)->K+K-)=1;  

  complex<double> ChiToA980pi0AmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToa0PiMag;
  
  for ( itChiToa0PiMag=ChiToa0PiMag.begin(); itChiToa0PiMag!=ChiToa0PiMag.end(); ++itChiToa0PiMag){
    boost::shared_ptr<const JPCLS> ChiToa0State=itChiToa0PiMag->first;
    double theChiToa0Mag=itChiToa0PiMag->second;
    double theChiToa0Phi=ChiToa0PiPhi[ChiToa0State];
    complex<double> expiphiChiToa0(cos(theChiToa0Phi), sin(theChiToa0Phi));
    
    ChiToA980pi0AmpTmpls+=theChiToa0Mag*expiphiChiToa0*sqrt(2*ChiToa0State->L+1); 
    //Clebsch(ChiToa0State->L=1,0,ChiToa0State->S=0, a0lam-pilam=0, ChiToa0State->J=1, 0)=1  
    //Clebsch(a0.J=0, a0lam=0, piJ=0, pilam=0, ChiToa0State->S=0, a0lam-pilam=0)=1
  }
  
  result*=ChiToA980pi0AmpTmpls;

  return result;
}

// complex<double> AbsPsi2SToKpKmPiGamLh::a2Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa2Pi, double a2Mass, double a2Width, Spin& lamChi){

//   complex<double> result(0.,0.);
//   Vector4<double> KpKm4V(theData->KpKm_HeliChic1_4V.E(), theData->KpKm_HeliChic1_4V.Px(), theData->KpKm_HeliChic1_4V.Py(), theData->KpKm_HeliChic1_4V.Pz());

//   for (Spin lam_a2=-2; lam_a2<=2; lam_a2++){
    
//     complex<double> ChiToa2PiAmpTmp1=conj(theData->DfChiToa2Pi[1][lamChi][lam_a2])
//       *sqrt(5)*conj(theData->Dfa2[2][lam_a2][0])*BreitWigner(KpKm4V, a2Mass, a2Width);
    
//     complex<double> ChiToa2PiAmpTmpls(0.,0.);
    
//     std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToa2Pi;
//     for ( itChiToa2Pi=ChiToa2Pi.begin(); itChiToa2Pi!=ChiToa2Pi.end(); ++itChiToa2Pi){
//          boost::shared_ptr<const JPCLS> ChiToa2PiState=itChiToa2Pi->first;
            
//          if ( fabs(lam_a2)>ChiToa2PiState->J || fabs(lam_a2)>ChiToa2PiState->S) continue;
            
//          double theChiToa2PiMag=itChiToa2Pi->second.first;
//          double theChiToa2PiPhi=itChiToa2Pi->second.second;
//          complex<double> expiphiChiToa2Pi(cos(theChiToa2PiPhi), sin(theChiToa2PiPhi));
            
//          ChiToa2PiAmpTmpls+=theChiToa2PiMag*expiphiChiToa2Pi*sqrt(2*ChiToa2PiState->L+1)
//            *Clebsch(ChiToa2PiState->L,0,ChiToa2PiState->S, lam_a2, ChiToa2PiState->J, lam_a2)
//            *Clebsch(2, lam_a2, 0, 0, ChiToa2PiState->S, lam_a2);
//           }
          
//        ChiToa2PiAmpTmp1*=ChiToa2PiAmpTmpls;
//        result+=ChiToa2PiAmpTmp1;
//   }

//   return result;
// }


complex<double> Psi2SToKpKmPiGamBaseLh::K0_Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK_0_KPhi , double K0_Mass, double K0_Width, Spin& lamChi){

  Vector4<double>& KpPi4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];

  Vector4<double>& KmPi4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];

  complex<double> result=conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKpPi][1][lamChi][0])
    *BreitWigner(KpPi4V, K0_Mass, K0_Width)
    +conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKmPi][1][lamChi][0])
    *BreitWigner(KmPi4V, K0_Mass, K0_Width);


  complex<double> ChiToK1400_0_KAmpTmpls(0.,0.);
  
  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToK1400_0_KMag;
  for ( itChiToK1400_0_KMag=ChiToK_0_KMag.begin(); itChiToK1400_0_KMag!=ChiToK_0_KMag.end(); ++itChiToK1400_0_KMag){
 
    boost::shared_ptr<const JPCLS> ChiToK1400_0_KState=itChiToK1400_0_KMag->first;
    double theChiToK1400_0_KMag=itChiToK1400_0_KMag->second;
    double theChiToK1400_0_KPhi=ChiToK_0_KPhi[ChiToK1400_0_KState];
    complex<double> expiphiChiToK1400_0_K(cos(theChiToK1400_0_KPhi), sin(theChiToK1400_0_KPhi));
    
    ChiToK1400_0_KAmpTmpls+=theChiToK1400_0_KMag*expiphiChiToK1400_0_K*sqrt(2*ChiToK1400_0_KState->L+1);
    //      *Clebsch(ChiToK1400_0_KState->L,0,ChiToK1400_0_KState->S, 0, ChiToK1400_0_KState->J, 0)
    //      *Clebsch(0, 0, 0, 0, 0, 0)
  }

  result*=ChiToK1400_0_KAmpTmpls;


  return result;
}


complex<double> Psi2SToKpKmPiGamBaseLh::K1_Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1_KPhi , double K1_Mass, double K1_Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  Vector4<double> KmPi4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];

  for (Spin lamK1=-1; lamK1<=1; lamK1++){
    
//     complex<double> ChiToK1_1400KAmpTmp1=conj(theData->DfChiToKpPi[1][lamChi][lamK1400])
//       *sqrt(3)*conj(theData->DfKp1[1][lamK1400][0])*BreitWigner(KpPi4V, K1_1400Mass, K1_1400Width)
//       +conj(theData->DfChiToKmPi[1][lamChi][lamK1400])
//       *sqrt(3)*conj(theData->DfKm1[1][lamK1400][0])*BreitWigner(KmPi4V, K1_1400Mass, K1_1400Width);

    complex<double> ChiToK1_KAmpTmp1=conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKpPi][1][lamChi][lamK1])
      *sqrt(3)*conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfKp1][1][lamK1][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K1_Mass, K1_Width, 1)
      +conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKmPi][1][lamChi][lamK1])
      *sqrt(3)*conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfKm1][1][lamK1][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K1_Mass, K1_Width, 1);
    complex<double> ChiToK1_KAmpTmpls(0.,0.);

    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToK1_KMag;
    for ( itChiToK1_KMag=ChiToK1_KMag.begin(); itChiToK1_KMag!=ChiToK1_KMag.end(); ++itChiToK1_KMag){
      boost::shared_ptr<const JPCLS> ChiToK1_KState=itChiToK1_KMag->first;
      double theChiToK1_KMag=itChiToK1_KMag->second;
      double theChiToK1_KPhi=ChiToK1_KPhi[ChiToK1_KState];
      complex<double> expiphiChiToK1_K(cos(theChiToK1_KPhi), sin(theChiToK1_KPhi));
      
      ChiToK1_KAmpTmpls+=theChiToK1_KMag*expiphiChiToK1_K*sqrt(2*ChiToK1_KState->L+1)
        *Clebsch(ChiToK1_KState->L,0,ChiToK1_KState->S, lamK1, ChiToK1_KState->J, lamK1)
        *Clebsch(1, lamK1, 0, 0, ChiToK1_KState->S, lamK1);
    }
    
    ChiToK1_KAmpTmp1*=ChiToK1_KAmpTmpls;
    result+=ChiToK1_KAmpTmp1;
  }

  return result;
}


complex<double> Psi2SToKpKmPiGamBaseLh::K2_Amp(EvtData* theData, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KMag, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToK1400_2_KPhi , double K2_1400Mass, double K2_1400Width, Spin& lamChi){

  complex<double> result(0.,0.);
  Vector4<double> KpPi4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KpPi_HeliChic1];
  Vector4<double> KmPi4V=theData->FourVecs[enumChic1ToKpKmPiGamData::KmPi_HeliChic1];

  for (Spin lamK1400_2=-2; lamK1400_2<=2; lamK1400_2++){
    
    complex<double> ChiToK1400_2_KAmpTmp1=conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKpPi][1][lamChi][lamK1400_2])
      *sqrt(5)*conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfKp2][2][lamK1400_2][0])*BreitWignerBlattW(KpPi4V, 0.493677, 0.1349766, K2_1400Mass, K2_1400Width, 2)
      +conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfChiToKmPi][1][lamChi][lamK1400_2])
      *sqrt(5)*conj(theData->WignerDs[enumChic1ToKpKmPiGamData::DfKm2][2][lamK1400_2][0])*BreitWignerBlattW(KmPi4V, 0.493677, 0.1349766, K2_1400Mass, K2_1400Width, 2);
    
    complex<double> ChiToK1400_2_KAmpTmpls(0.,0.);
    
    std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >::iterator itChiToK1400_2_KMag;
    for ( itChiToK1400_2_KMag=ChiToK1400_2_KMag.begin(); itChiToK1400_2_KMag!=ChiToK1400_2_KMag.end(); ++itChiToK1400_2_KMag){
            boost::shared_ptr<const JPCLS> ChiToK1400_2_KState=itChiToK1400_2_KMag->first;
            
            if ( fabs(lamK1400_2)>ChiToK1400_2_KState->J || fabs(lamK1400_2)>ChiToK1400_2_KState->S) continue;
            
            double theChiToK1400_2_KMag=itChiToK1400_2_KMag->second;
            double theChiToK1400_2_KPhi=ChiToK1400_2_KPhi[ChiToK1400_2_KState];
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



// complex<double> AbsPsi2SToKpKmPiGamLh::KKPi_Amp(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi, Spin& lamChi){

//   complex<double> result=conj(theData->DfChiToKKPi[1][lamChi][0]);   

//   complex<double> ChiToKKPiAmpTmpls(0.,0.);
  
//   std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itChiToKKPi;
  
//   for ( itChiToKKPi=ChiToKKPi.begin(); itChiToKKPi!=ChiToKKPi.end(); ++itChiToKKPi){
//     boost::shared_ptr<const JPCLS> ChiToKKPiState=itChiToKKPi->first;
//     double theChiToKKPiMag=itChiToKKPi->second.first;
//     double theChiToKKPiPhi=itChiToKKPi->second.second;
//     complex<double> expiphiChiToKKPi(cos(theChiToKKPiPhi), sin(theChiToKKPiPhi));
    
//     ChiToKKPiAmpTmpls+=theChiToKKPiMag*expiphiChiToKKPi*sqrt(2.*ChiToKKPiState->L+1); 
//     //Clebsch(ChiToKKPiState->L=1,0,ChiToKKPiState->S=0, lam=0, ChiToKKPiState->J=1, 0)=1  
//     //Clebsch(0, 0, 0, 0, ChiToKKPiState->S=0, 0)=1
//   }
  
//   result*=ChiToKKPiAmpTmpls;

//   return result;
// }
