#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKBaseLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKEventList.hh"
#include "ErrLogger/ErrLogger.hh"

JpsiGamKsKlKKBaseLh::JpsiGamKsKlKKBaseLh(boost::shared_ptr<const JpsiGamKsKlKKEventList> theEvtList, const std::map<const std::string, bool>& hypMap) :
  AbsJpsiGamKsKlKKLh(theEvtList)
  ,_K0_1430Hyp(true)
  ,_K1_1410Hyp(true)
  ,_K2_1430Hyp(true)
  ,_K1_1680Hyp(false)
  ,_KKPi_Hyp(false)
  ,_nFitParams(0)  
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
  _ampVec.push_back(paramEnumKpKmPiGam::ChiGam);
  _ampVec.push_back(paramEnumKpKmPiGam::K890K);
  if(_K0_1430Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_0_1400K);
  if(_K1_1410Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_1_1400K);
  if(_K2_1430Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_2_1400K);
  if(_K1_1680Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_1_1680K);
  if(_KKPi_Hyp) _ampVec.push_back(paramEnumKpKmPiGam::KKPi);
  _ampVec.push_back(paramEnumKpKmPiGam::a980Pi);

  _massVec.push_back(paramEnumKpKmPiGam::K890);
  if(_K0_1430Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_0_1400);
  if(_K1_1410Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_1_1400);
  if(_K2_1430Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_2_1400);
  if(_K1_1680Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_1_1680);
  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  std::vector<unsigned int>::iterator massIt; 
  for (massIt=_massVec.begin(); massIt!=_massVec.end(); ++massIt){
    _nFitParams+=2;
  }

  _nFitParams+=3; //a(980) Flatte parameters
  _nFitParams+=1; //phase space
}

JpsiGamKsKlKKBaseLh::JpsiGamKsKlKKBaseLh( boost::shared_ptr<AbsJpsiGamKsKlKKLh> theLhPtr, const std::map<const std::string, bool>& hypMap ) :
  AbsJpsiGamKsKlKKLh(theLhPtr->getEventList())
  ,_K0_1430Hyp(true)
  ,_K1_1410Hyp(true)
  ,_K2_1430Hyp(true)
  ,_K1_1680Hyp(false)
  ,_KKPi_Hyp(false)
  ,_nFitParams(0)
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

  _ampVec.push_back(paramEnumKpKmPiGam::ChiGam);
  _ampVec.push_back(paramEnumKpKmPiGam::K890K);
  if(_K0_1430Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_0_1400K);
  if(_K1_1410Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_1_1400K);
  if(_K2_1430Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_2_1400K);
  if(_K1_1680Hyp) _ampVec.push_back(paramEnumKpKmPiGam::K_1_1680K);
  if(_KKPi_Hyp) _ampVec.push_back(paramEnumKpKmPiGam::KKPi);
  _ampVec.push_back(paramEnumKpKmPiGam::a980Pi);

  _massVec.push_back(paramEnumKpKmPiGam::K890);
  if(_K0_1430Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_0_1400);
  if(_K1_1410Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_1_1400);
  if(_K2_1430Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_2_1400);
  if(_K1_1680Hyp) _massVec.push_back(paramEnumKpKmPiGam::K_1_1680);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }

  std::vector<unsigned int>::iterator massIt; 
  for (massIt=_massVec.begin(); massIt!=_massVec.end(); ++massIt){
    _nFitParams+=2;
  }

  _nFitParams+=3; //a(980) Flatte parameters
  _nFitParams+=1; //phase space parameter
}

JpsiGamKsKlKKBaseLh::~JpsiGamKsKlKKBaseLh()
{;
}


complex<double> JpsiGamKsKlKKBaseLh::calcCoherentAmp(Spin Minit, Spin lamGam, const paramKpKmPiGam& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){

  complex<double> result(0.,0.);

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;

  for (Spin lamChi=-1; lamChi<=1; lamChi++){

      Spin lamChiGam=lamChi-lamGam;

      complex<double> PsiDecAmp=conj(theData->dfPsi[1][Minit][lamChiGam]);
      complex<double> PsiDecAmpTmpls(0.,0.);  
      
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
      
      for ( itPsi=PsiToChiGam.begin(); itPsi!=PsiToChiGam.end(); ++itPsi){
    
	boost::shared_ptr<const JPCLS> PsiState=itPsi->first;

        if (fabs(lamChiGam)>PsiState->J || fabs(lamChiGam)>PsiState->S) continue;
	
	double thePsiMag=itPsi->second.first;
	double thePsiPhi=itPsi->second.second;
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

complex<double> JpsiGamKsKlKKBaseLh::calcDecAmp(complex<double>& inAmp,Spin lamChi, const paramKpKmPiGam& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData){

  complex<double> result=inAmp;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK890K=theParamVal.ChiToK890K;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToa0Pi=theParamVal.ChiToa0Pi;
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi=theParamVal.ChiToKKPi;

  double K890Mass=theParamVal.BwK890.first;
  double K890Width=theParamVal.BwK890.second;
  double a980FlatteMass=theParamVal.FlatMa980;
  double a980FlatgKK=theParamVal.FlatgKK;
  double a980FlatgEtaPi=theParamVal.FlatgEtaPi;

  //Chi to a0(980) pi0 and a0(980)->K+K-
  result+=a980FlatteAmp(theData, ChiToa0Pi, a980FlatteMass, a980FlatgKK, a980FlatgEtaPi, lamChi);
  
  //Chi -> K890 K and K890->K pi0
  result+=K892Amp(theData, ChiToK890K, K890Mass, K890Width, lamChi);
  
  //Chi -> K_0_1400 K and K_0_1400->K pi0
  if(_K0_1430Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_0_1400K=theParamVal.ChiToK1400_0_K;
    double K_0_1400Mass=theParamVal.BwK1400_0.first;
    double K_0_1400Width=theParamVal.BwK1400_0.second;
    result+=K0_1400Amp(theData, ChiToK_0_1400K, K_0_1400Mass, K_0_1400Width, lamChi);
  }
  
  //Chi -> K_1_1400 K and K_1_1400->K pi0
  if(_K1_1410Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1400K=theParamVal.ChiToK1400_1_K; 
    double K_1_1400Mass=theParamVal.BwK1400_1.first;
    double K_1_1400Width=theParamVal.BwK1400_1.second;
    result+=K1_1400Amp(theData, ChiToK_1_1400K, K_1_1400Mass, K_1_1400Width, lamChi);
  }
  
  //Chi -> K_2_1400 K and K_2_1400->K pi0
  if(_K2_1430Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_2_1400K=theParamVal.ChiToK1400_2_K;
    double K_2_1400Mass=theParamVal.BwK1400_2.first;
    double K_2_1400Width=theParamVal.BwK1400_2.second;
    result+=K2_1400Amp(theData, ChiToK_2_1400K, K_2_1400Mass, K_2_1400Width, lamChi);
  }


  //K_1_1680->K pi0
  if(_K1_1680Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToK_1_1680K=theParamVal.ChiToK1680_1_K; 
    double K_1_1680Mass=theParamVal.BwK1680_1.first;
    double K_1_1680Width=theParamVal.BwK1680_1.second;
    result+=K892Amp(theData, ChiToK_1_1680K, K_1_1680Mass, K_1_1680Width, lamChi);
  }

  //Chi -> K+ K- pi0
  if(_KKPi_Hyp){
    std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ChiToKKPi=theParamVal.ChiToKKPi;
    result+=KKPi_Amp(theData, ChiToKKPi, lamChi);
  }
  
  return result;
}


void JpsiGamKsKlKKBaseLh::setMnUsrParams(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal){

  _fitParamsKpKmPiGam.setMnUsrParamsFlattea980Mass(upar, startVal, errVal, paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::phaseSpace));
  
  _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal,paramEnumKpKmPiGam::ChiGam);
  _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::K890K);

  if(_K0_1430Hyp) _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::K_0_1400K);
  if(_K1_1410Hyp) _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::K_1_1400K);
  if(_K2_1430Hyp) _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::K_2_1400K);
  if(_K1_1680Hyp) _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::K_1_1680K);
  if(_KKPi_Hyp) _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::KKPi); 
 
  _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnumKpKmPiGam::a980Pi);

//K*(982) BW
  _fitParamsKpKmPiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::K890);
   
//K*(1400) BWs
  if(_K0_1430Hyp) _fitParamsKpKmPiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::K_0_1400);  
  if(_K1_1410Hyp) _fitParamsKpKmPiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::K_1_1400);   
  if(_K2_1430Hyp) _fitParamsKpKmPiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::K_2_1400);
  if(_K1_1680Hyp) _fitParamsKpKmPiGam.setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::K_1_1680);   
     
//a(980) BW
//    setMnUsrParamsMass(upar, startVal, errVal, paramEnumKpKmPiGam::a980"); 
   _fitParamsKpKmPiGam.setMnUsrParamsFlattea980Mass(upar, startVal, errVal, paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::a980));
   
}


int JpsiGamKsKlKKBaseLh::setFitParamVal(paramKpKmPiGam& theParamVal, const std::vector<double>& par) {
  
  if (par.size()!= nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  int counter=0;

  counter=_fitParamsKpKmPiGam.setFitParamFlattea980Mass(theParamVal, par, counter, paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::phaseSpace));

  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::ChiGam);

  //K*(982) amplitude params
  counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::K890K);

  //K*(1400)_1 amplitude params
  if(_K0_1430Hyp) counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::K_0_1400K);
  if(_K1_1410Hyp) counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::K_1_1400K);
  if(_K2_1430Hyp) counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::K_2_1400K);
  if(_K1_1680Hyp) counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::K_1_1680K);
  if(_KKPi_Hyp) counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::KKPi);

  //a(980) amplitude params
  counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::a980Pi);

  // K890 mass 
  counter=_fitParamsKpKmPiGam.setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::K890);

  // K1400 masses
  if(_K0_1430Hyp) counter=_fitParamsKpKmPiGam.setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::K_0_1400);
  if(_K1_1410Hyp) counter=_fitParamsKpKmPiGam.setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::K_1_1400);
  if(_K2_1430Hyp) counter=_fitParamsKpKmPiGam.setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::K_2_1400);
  if(_K1_1680Hyp) counter=_fitParamsKpKmPiGam.setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::K_1_1680);
  // a(980) mass 
//   counter=setFitParamValMass(theParamVal, par, counter, paramEnumKpKmPiGam::a980);
  counter=_fitParamsKpKmPiGam.setFitParamFlattea980Mass(theParamVal, par, counter, paramEnumKpKmPiGam::name(paramEnumKpKmPiGam::a980));
  return counter;
}

unsigned int  JpsiGamKsKlKKBaseLh::nFitParams(){

  return _nFitParams;
}


void JpsiGamKsKlKKBaseLh::print(std::ostream& os) const{
  os << "JpsiGamKsKlKKBaseLh::print\n";
}

void JpsiGamKsKlKKBaseLh::printCurrentFitResult(paramKpKmPiGam& theParamVal) {

  DebugMsg<< "phaseSpace: " << theParamVal.phaseSpace << endmsg;


  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< paramEnumKpKmPiGam::name(*itAmps) << endmsg;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParamsKpKmPiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }  
  }


  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    DebugMsg<< paramEnumKpKmPiGam::name(*itMasses) << endmsg;
    std::pair<double, double> tmpParam=_fitParamsKpKmPiGam.massPair(theParamVal, *itMasses);
    DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
  }

  DebugMsg<< "a980:" << endmsg;
  DebugMsg <<"\t mass:" << theParamVal.FlatMa980 <<"\t gKK:" <<  theParamVal.FlatgKK  << "\t gEtaPi:" <<  theParamVal.FlatgEtaPi << endmsg;
}

void JpsiGamKsKlKKBaseLh::dumpCurrentResult(std::ostream& os, paramKpKmPiGam& theParamVal, std::string& suffix){

  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  std::string tmpStringPs="phaseSpace"+suffix;
  os << tmpStringPs << "\t" << theParamVal.phaseSpace  << "\t" << 0. << std::endl;


  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParamsKpKmPiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];

      std::string tmpStringDec=(*itJPCLS)->name()+paramEnumKpKmPiGam::name(*itAmps)+suffix;
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }  
  }

  std::vector<unsigned int>::const_iterator itMasses;
  for ( itMasses=_massVec.begin(); itMasses!=_massVec.end(); ++itMasses){
    std::string tmpStringMass=paramEnumKpKmPiGam::name(*itMasses)+suffix;

    std::pair<double, double> tmpParam=_fitParamsKpKmPiGam.massPair(theParamVal, *itMasses);
    os << tmpStringMass << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
  }

  std::string tmpStringFlat="FlatMa980"+suffix;
  os << tmpStringFlat << "\t" << theParamVal.FlatMa980  << "\t" << 0. << std::endl;

  tmpStringFlat="FlatgKK"+suffix;
  os << tmpStringFlat << "\t" << theParamVal.FlatgKK  << "\t" << 0. << std::endl;

  tmpStringFlat="FlatgEtaPi"+suffix;
  os << tmpStringFlat << "\t" << theParamVal.FlatgEtaPi  << "\t" << 0. << std::endl;
}
