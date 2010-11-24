#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/Psi2STo2K2PiGam/PsiToChic0GamProdLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamEvtList.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "ErrLogger/ErrLogger.hh"

PsiToChic0GamProdLh::PsiToChic0GamProdLh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList> theEvtList, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates) :
  AbsPsi2STo2K2PiGamLh(theEvtList,theStates)  
{
}

PsiToChic0GamProdLh::PsiToChic0GamProdLh( boost::shared_ptr<AbsPsi2STo2K2PiGamLh> theLhPtr ) :
  AbsPsi2STo2K2PiGamLh(theLhPtr->getEventList(), theLhPtr->getPsi2STo2K2PiGamStates())
{
}

PsiToChic0GamProdLh::~PsiToChic0GamProdLh()
{;
}

double PsiToChic0GamProdLh::calcEvtIntensity(Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData, const Psi2STo2K2PiGamData::fitParamVal& theParamVal){

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

complex<double> PsiToChic0GamProdLh::chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  complex<double> result(0.,0.);
  return result;  //dummy !!!
}


complex<double> PsiToChic0GamProdLh::calcCoherentAmp(Spin Minit, Spin lamGam, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  complex<double> dummyresult(0.,0.);
  return dummyresult; 
}

complex<double> PsiToChic0GamProdLh::calcCoherentAmp(Spin Minit, Spin MChi, Spin MGamma, const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData){
  
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



void PsiToChic0GamProdLh::setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal){
  checkFitParamVal(startVal);
  checkFitParamVal(errVal);

  setMnUsrParamsDec(upar,startVal,errVal,"ChiGam");
}



int PsiToChic0GamProdLh::setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const{
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_Psi2STo2K2PiGamStatesPtr->PsiToChiGamStates();
  
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

void PsiToChic0GamProdLh::print(std::ostream& os) const{
  os << "PsiToChic0GamProdLh::print\n";
}


