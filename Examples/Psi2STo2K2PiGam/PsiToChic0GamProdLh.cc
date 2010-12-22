#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/PsiToChic0GamProdLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToChic0GamProdLh::PsiToChic0GamProdLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList) :
  AbsPsi2STo2K2PiGamLh(theEvtList)
  ,_nFitParams(0)  
{
  _ampVec.push_back(paramEnum2K2PiGam::ChiGam);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }
}

PsiToChic0GamProdLh::PsiToChic0GamProdLh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr ) :
  AbsPsi2STo2K2PiGamLh(theLhPtr->getEventList())
  ,_nFitParams(0)  
{
  _ampVec.push_back(paramEnum2K2PiGam::ChiGam);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }
}

PsiToChic0GamProdLh::~PsiToChic0GamProdLh()
{
  ;
}

double PsiToChic0GamProdLh::calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const param2K2PiGam& theParamVal){

  double result=0.;
  
   for (Spin PsiM=-1; PsiM<=1; PsiM=PsiM+2){
    for (Spin ChiM=-1; ChiM<=1; ChiM++){
      for (Spin GamM=-1; GamM<=1; GamM=GamM+2){

	result+=norm(calcCoherentAmp(PsiM, ChiM, GamM, theParamVal, theData));
      }
    }
  }
  return result;  
}

complex<double> PsiToChic0GamProdLh::chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  complex<double> result(0.,0.);
  return result;  //dummy !!!
}


complex<double> PsiToChic0GamProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}

complex<double> PsiToChic0GamProdLh::calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  
  complex<double> result(0.,0.);
  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > PsiToChiGam=theParamVal.PsiToChiGam;

  std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >::iterator itPsi;
  for ( itPsi=PsiToChiGam.begin(); itPsi!=PsiToChiGam.end(); ++itPsi){
    
    boost::shared_ptr<const JPCLS> PsiState=itPsi->first;
    double thePsiMag=itPsi->second.first;
    double thePsiPhi=itPsi->second.second;
    complex<double> expiphiPsi(cos(thePsiPhi), sin(thePsiPhi));
    
    Spin lamChiGam=MChi-MGamma;
    if (fabs(lamChiGam)>PsiState->J || fabs(lamChiGam)>PsiState->S) continue;
    
    complex<double> PsiDecAmp=thePsiMag*expiphiPsi*sqrt(2*PsiState->L+1)
      *Clebsch(PsiState->L,0,PsiState->S,lamChiGam,PsiState->J,lamChiGam)
      *Clebsch(1, MChi, 1, -MGamma, PsiState->S,lamChiGam)
      *conj(theData->DfPsi[PsiState->J][Minit][lamChiGam]);
    
    result+=PsiDecAmp;
  }      

  return result; 
}



void PsiToChic0GamProdLh::setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal, param2K2PiGam& errVal){
  _fitParams2K2PiGam.setMnUsrParamsDec(upar,startVal,errVal, paramEnum2K2PiGam::ChiGam);
}



int PsiToChic0GamProdLh::setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par){

  if (par.size() != nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }  

  int counter=0;

  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter= _fitParams2K2PiGam.setFitParamValDec(theParamVal, par, counter, paramEnum2K2PiGam::ChiGam);

  return counter;
}

void PsiToChic0GamProdLh::print(std::ostream& os) const{
  os << "PsiToChic0GamProdLh::print\n";
}

void PsiToChic0GamProdLh::printCurrentFitResult(param2K2PiGam& theParamVal){

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;  
  
  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      DebugMsg<< (*itJPCLS)->name()<< paramEnum2K2PiGam::name(*itAmps) << endmsg;
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParams2K2PiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      DebugMsg <<"\t mag:" << tmpParam.first <<"\t phi:" << tmpParam.second  << endmsg;
    }  
  }
}

void PsiToChic0GamProdLh::dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix){
  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;

  std::vector<unsigned int>::const_iterator itAmps;
  for ( itAmps=_ampVec.begin(); itAmps!=_ampVec.end(); ++itAmps){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParams2K2PiGam.jpclsVec(*itAmps);
    
    for ( itJPCLS=JPCLSs.begin(); itJPCLS!=JPCLSs.end(); ++itJPCLS){
      std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > currentMap=_fitParams2K2PiGam.ampMap(theParamVal, *itAmps);
      std::pair<double, double> tmpParam=currentMap[(*itJPCLS)];
      
      std::string tmpStringDec=(*itJPCLS)->name()+paramEnum2K2PiGam::name(*itAmps)+suffix;
      os << tmpStringDec << "\t" << tmpParam.first  << "\t" << tmpParam.second << std::endl;
    }  
  }
}
