#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamBaseLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "ErrLogger/ErrLogger.hh"

Psi2SToKpKmPiGamBaseLh::Psi2SToKpKmPiGamBaseLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsPsi2SToKpKmPiGamLh(theEvtList)
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

Psi2SToKpKmPiGamBaseLh::Psi2SToKpKmPiGamBaseLh( boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsPsi2SToKpKmPiGamLh(theLhPtr->getEventList())
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


complex<double> Psi2SToKpKmPiGamBaseLh::calcCoherentAmp(Spin Minit, Spin lamGam, fitParams& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){

  complex<double> result(0.,0.);

  //  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& thePsiMags=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::ChiGam];
      std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& thePsiPhis=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::ChiGam]; 

  for (Spin lamChi=-1; lamChi<=1; lamChi++){

      Spin lamChiGam=lamChi-lamGam;

      complex<double> PsiDecAmp=conj(theData->dfPsi[1][Minit][lamChiGam]);
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

complex<double> Psi2SToKpKmPiGamBaseLh::calcDecAmp(complex<double>& inAmp,Spin lamChi, fitParams& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){

  complex<double> result=inAmp;
//   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToK890KMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::K890K];
//   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToK890KPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::K890K];
//   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::a980Pi];
//   std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > ChiToa0PiPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::a980Pi];
//  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToKKPiMag=theParamVal.Mags[paramEnumChic1ToKpKmPiGam::KKPi];
//  std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& ChiToKKPiPhi=theParamVal.Phis[paramEnumChic1ToKpKmPiGam::KKPi];
//   double K890Mass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::K890];
//   double K890Width=theParamVal.Widths[paramEnumChic1ToKpKmPiGam::K890];
//   double a980FlatteMass=theParamVal.Masses[paramEnumChic1ToKpKmPiGam::a980];
//   double a980FlatgKK=theParamVal.gFactors[paramEnumChic1ToKpKmPiGam::a980gKK];
//   double a980FlatgEtaPi=theParamVal.gFactors[paramEnumChic1ToKpKmPiGam::a980gEtaPi];

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

