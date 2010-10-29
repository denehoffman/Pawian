#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/K1400Spin01Lh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

K1400Spin01Lh::K1400Spin01Lh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList, boost::shared_ptr<const Psi2SToKpKmPiGamStates> theStates) :
  AbsPsi2SToKpKmPiGamLh(theEvtList,theStates)  
{
}

K1400Spin01Lh::K1400Spin01Lh( boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theLhPtr ) :
  AbsPsi2SToKpKmPiGamLh(theLhPtr->getEventList(), theLhPtr->getPsi2SToKpKmPiGamStates())
{
}

K1400Spin01Lh::~K1400Spin01Lh()
{;
}

complex<double> K1400Spin01Lh::calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){

  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K=theParamVal.ChiToK890K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_0_K=theParamVal.ChiToK1400_0_K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK1400_1_K=theParamVal.ChiToK1400_1_K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi=theParamVal.ChiToa0Pi;

  double K890Mass=theParamVal.BwK890.first;
  double K890Width=theParamVal.BwK890.second;
  double K1400_0_Mass=theParamVal.BwK1400_0.first;
  double K1400_0_Width=theParamVal.BwK1400_0.second;
  double K1400_1_Mass=theParamVal.BwK1400_1.first;
  double K1400_1_Width=theParamVal.BwK1400_1.second;
//   double a980Mass=theParamVal.Bwa980.first;
//   double a980Width=theParamVal.Bwa980.second;
  double a980FlatteMass=theParamVal.FlatMa980;
  double a980FlatgKK=theParamVal.FlatgKK;
  double a980FlatgEtaPi=theParamVal.FlatgEtaPi;

  for (Spin lamChi=-1; lamChi<=1; lamChi++){
      
      Spin lamChiGam=lamChi-lamGam;
      complex<double> PsiDecAmp=conj(theData->DfPsi[1][Minit][lamChiGam]);
      complex<double> PsiDecAmpTmpls(0.,0.);


      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
      for ( itPsi=PsiToChiGam.begin(); itPsi!=PsiToChiGam.end(); ++itPsi){
	
	boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
	
        if (fabs(lamChiGam)>PsiState->J || fabs(lamChiGam)>PsiState->S) continue;
	
	double thePsiMag=itPsi->second.first;
	double thePsiPhi=itPsi->second.second;
	complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
    

	//Chi to a0(980) pi0 and a0(980)->K+K-
// 	complex<double> ChiToA980pi0Amp=a980Amp(theData, ChiToa0Pi, a980Mass, a980Width, lamChi);
	complex<double> ChiToA980pi0Amp=a980FlatteAmp(theData, ChiToa0Pi, a980FlatteMass, a980FlatgKK, a980FlatgEtaPi, lamChi);

	//Chi -> K890 K and K890->K pi0
	complex<double> ChiToK890KAmp=K892Amp(theData, ChiToK890K, K890Mass, K890Width, lamChi);


	//Chi -> K1400_0 K and K1400_0->K pi0
	complex<double> ChiToK1400_0_KAmp=K0_1400Amp(theData, ChiToK1400_0_K, K1400_0_Mass, K1400_0_Width, lamChi);

        //Chi -> K1400_1 K and K1400_1->K pi0
        complex<double> ChiToK1400_1_KAmp=K1_1400Amp(theData, ChiToK1400_1_K, K1400_1_Mass, K1400_1_Width, lamChi);


	PsiDecAmpTmpls+=thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
          *Clebsch(PsiState->L,0,PsiState->S,lamChiGam,PsiState->J,lamChiGam)
          *Clebsch(1, lamChi, 1, -lamGam, PsiState->S,lamChiGam)
          *(ChiToA980pi0Amp+ChiToK890KAmp+ChiToK1400_0_KAmp+ChiToK1400_1_KAmp);
      }
      
      PsiDecAmp*=PsiDecAmpTmpls;
      result+=PsiDecAmp;
    
  }
      return result;
}


void K1400Spin01Lh::setMnUsrParams(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);
//   setMnUsrParamsPsiDec(upar, startVal, errVal);
  setMnUsrParamsDec(upar,startVal,errVal,"ChiGam");
  setMnUsrParamsDec(upar,startVal,errVal, "K890Pi");
  setMnUsrParamsDec(upar,startVal,errVal, "K_0_1400Pi");
  setMnUsrParamsDec(upar,startVal,errVal, "K_1_1400Pi");
  setMnUsrParamsDec(upar,startVal,errVal, "a980Pi");

//K*(982) BW
  setMnUsrParamsMass(upar, startVal, errVal, "K890");
   
//K*(1400)_0 BW 
  setMnUsrParamsMass(upar, startVal, errVal, "K_0_1400"); 

//K*(1400)_1 BW
  setMnUsrParamsMass(upar, startVal, errVal, "K_1_1400");   
   
//a(980) BW
//    setMnUsrParamsMass(upar, startVal, errVal, "a980");
  setMnUsrParamsFlattea980Mass(upar, startVal, errVal, "a980");
   
}


int K1400Spin01Lh::setFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2SToKpKmPiGamStatesPtr->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToKst0K=_Psi2SToKpKmPiGamStatesPtr->ChiToKst0KStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToKst1K=_Psi2SToKpKmPiGamStatesPtr->ChiToKst1KStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToa980pi=_Psi2SToKpKmPiGamStatesPtr->ChiToa0PiStates();
  
  if (par.size()!= JPCLSPsiToChiGam.size()*2-1+JPCLSChiToKst0K.size()*2+JPCLSChiToKst1K.size()*4+JPCLSChiToa980pi.size()*2+9) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << JPCLSPsiToChiGam.size()*2+1+JPCLSChiToKst0K.size()*2+JPCLSChiToKst1K.size()*4+JPCLSChiToa980pi.size()*2+9 << endmsg;
    exit(1);
  }

  int counter=0;

  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "ChiGam");

  //K*(982) amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "K890Pi");

  //K*(1400)_0 amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "K_0_1400Pi");

  //K*(1400)_1 amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "K_1_1400Pi");

  //a(980) amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "a980Pi");

  // K890 mass 
  counter=setFitParamValMass(theParamVal, par, counter, "K890");

  // K1400_0 mass 
  counter=setFitParamValMass(theParamVal, par, counter, "K_0_1400");

  // K1400_2 mass 
  counter=setFitParamValMass(theParamVal, par, counter, "K_1_1400");

  // a(980) mass 
//   counter=setFitParamValMass(theParamVal, par, counter, "a980");
  counter=setFitParamFlattea980Mass(theParamVal, par, counter, "Flattea980");
  return counter;
}

void K1400Spin01Lh::print(std::ostream& os) const{
  os << "K1400Spin01Lh::print\n";
}

void K1400Spin01Lh::printCurrentFitResult(Psi2SToKpKmPiGamData::fitParamVal& theParamVal) const{

  AbsPsi2SToKpKmPiGamLh::printCurrentFitResult(theParamVal);
  
  //  print fit paramss
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToKst0K=_Psi2SToKpKmPiGamStatesPtr->ChiToKst0KStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToKst1K=_Psi2SToKpKmPiGamStatesPtr->ChiToKst1KStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToa980pi=_Psi2SToKpKmPiGamStatesPtr->ChiToa0PiStates();


  for ( itJPCLS=JPCLSChiToKst1K.begin(); itJPCLS!=JPCLSChiToKst1K.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "K890Pi" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiToK890K[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  for ( itJPCLS=JPCLSChiToKst0K.begin(); itJPCLS!=JPCLSChiToKst0K.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "K1400_0_Pi" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiToK1400_0_K[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  for ( itJPCLS=JPCLSChiToKst1K.begin(); itJPCLS!=JPCLSChiToKst1K.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "K1400_1_Pi" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiToK1400_1_K[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  for ( itJPCLS=JPCLSChiToa980pi.begin(); itJPCLS!=JPCLSChiToa980pi.end(); ++itJPCLS){
    DebugMsg<< (*itJPCLS)->name()<< "a980Pi" << endmsg;
    std::pair<double, double> tmpParam=theParamVal.ChiToa0Pi[(*itJPCLS)];
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

    DebugMsg<< "K890:" << endmsg;
    std::pair<double, double> tmpParamK890=theParamVal.BwK890;
    DebugMsg <<"\t mass:" << tmpParamK890.first <<"\t width:" << tmpParamK890.second  << endmsg;

    DebugMsg<< "K1400_0:" << endmsg;
    std::pair<double, double> tmpParamK1400_0=theParamVal.BwK1400_0;
    DebugMsg <<"\t mass:" << tmpParamK1400_0.first <<"\t width:" << tmpParamK1400_0.second  << endmsg;

    DebugMsg<< "K1400_1:" << endmsg;
    std::pair<double, double> tmpParamK1400_1=theParamVal.BwK1400_1;
    DebugMsg <<"\t mass:" << tmpParamK1400_1.first <<"\t width:" << tmpParamK1400_1.second  << endmsg;

//     DebugMsg<< "a980:" << endmsg;
//     std::pair<double, double> tmpParama980=theParamVal.Bwa980;
//     DebugMsg <<"\t mass:" << tmpParama980.first <<"\t width:" << tmpParama980.second  << endmsg;
    DebugMsg<< "a980:" << endmsg;
    DebugMsg <<"\t mass:" << theParamVal.FlatMa980 <<"\t gKK:" <<  theParamVal.FlatgKK  << "\t gEtaPi:" <<  theParamVal.FlatgEtaPi << endmsg;
}
