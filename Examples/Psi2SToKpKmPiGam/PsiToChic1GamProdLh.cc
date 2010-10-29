#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2SToKpKmPiGam/PsiToChic1GamProdLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamEventList.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToChic1GamProdLh::PsiToChic1GamProdLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList> theEvtList, boost::shared_ptr<const Psi2SToKpKmPiGamStates> theStates) :
  AbsPsi2SToKpKmPiGamLh(theEvtList,theStates)  
{
}

PsiToChic1GamProdLh::PsiToChic1GamProdLh( boost::shared_ptr<AbsPsi2SToKpKmPiGamLh> theLhPtr ) :
  AbsPsi2SToKpKmPiGamLh(theLhPtr->getEventList(), theLhPtr->getPsi2SToKpKmPiGamStates())
{
}

PsiToChic1GamProdLh::~PsiToChic1GamProdLh()
{;
}

double PsiToChic1GamProdLh::calcEvtIntensity(Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal){

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

complex<double> PsiToChic1GamProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}

complex<double> PsiToChic1GamProdLh::calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const Psi2SToKpKmPiGamData::fitParamVal& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData){
  
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



void PsiToChic1GamProdLh::setMnUsrParams(MnUserParameters& upar, Psi2SToKpKmPiGamData::fitParamVal& startVal,  Psi2SToKpKmPiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  setMnUsrParamsDec(upar,startVal,errVal,"ChiGam");
}



int PsiToChic1GamProdLh::setFitParamVal(Psi2SToKpKmPiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2SToKpKmPiGamStatesPtr->PsiToChiGamStates();
  
  if (par.size()!= JPCLSPsiToChiGam.size()*2-1) {
    Alert << "size of parameters wrong!!! par.size()=" << par.size() << 
      "\t it should be" << JPCLSPsiToChiGam.size()*2-1 << endmsg;
    exit(1);
  }

  int counter=0;

  //Psi(2S) ->Chi_c1 gamma amplitude params
  counter=setFitParamValDec(theParamVal, par, counter, "ChiGam");

  return counter;
}

void PsiToChic1GamProdLh::print(std::ostream& os) const{
  os << "PsiToChic1GamProdLh::print\n";
}


