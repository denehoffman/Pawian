#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/PsiToChic1GamProdLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToChic1GamProdLh::PsiToChic1GamProdLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList) :
  AbsPsi2SToKpKmPiGamLh(theEvtList)
  ,_nFitParams(0)    
{
  _ampVec.push_back(paramEnumKpKmPiGam::ChiGam);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }
}

PsiToChic1GamProdLh::PsiToChic1GamProdLh( boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theLhPtr ) :
  AbsPsi2SToKpKmPiGamLh(theLhPtr->getEventList())
  ,_nFitParams(0) 
{
  _ampVec.push_back(paramEnumKpKmPiGam::ChiGam);

  std::vector<unsigned int>::iterator ampIt;
  for (ampIt=_ampVec.begin(); ampIt!=_ampVec.end(); ++ampIt){
    std::vector< boost::shared_ptr<const JPCLS> > JPCLSs=_fitParamsKpKmPiGam.jpclsVec(*ampIt);
    _nFitParams+=2*JPCLSs.size();
  }
}

PsiToChic1GamProdLh::~PsiToChic1GamProdLh()
{;
}

double PsiToChic1GamProdLh::calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, const paramKpKmPiGam& theParamVal){

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

complex<double> PsiToChic1GamProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}

complex<double> PsiToChic1GamProdLh::calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){
  
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



void PsiToChic1GamProdLh::setMnUsrParams(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal){

  _fitParamsKpKmPiGam.setMnUsrParamsDec(upar,startVal,errVal,paramEnumKpKmPiGam::ChiGam);
}



int PsiToChic1GamProdLh::setFitParamVal(paramKpKmPiGam& theParamVal, const std::vector<double>& par){
  
  if (par.size()!= nFitParams() ) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << nFitParams() << endmsg;
    exit(1);
  }

  int counter=0;

  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=_fitParamsKpKmPiGam.setFitParamValDec(theParamVal, par, counter, paramEnumKpKmPiGam::ChiGam);

  return counter;
}

unsigned int  PsiToChic1GamProdLh::nFitParams(){
  return _nFitParams;
}

void PsiToChic1GamProdLh::print(std::ostream& os) const{
  os << "PsiToChic1GamProdLh::print\n";
}

void PsiToChic1GamProdLh::printCurrentFitResult(paramKpKmPiGam& theParamVal){

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

}

void PsiToChic1GamProdLh::dumpCurrentResult(std::ostream& os, paramKpKmPiGam& theParamVal, std::string& suffix){

  if ( suffix.compare("Val") != 0 && suffix.compare("Err") !=0 ){
    Warning << "suffix " << suffix << " not supported!!! Use Val or Err" << endmsg;
    return;
  }

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

}
