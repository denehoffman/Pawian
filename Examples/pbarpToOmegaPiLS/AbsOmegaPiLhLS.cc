#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiEventListLS.hh"
#include "PwaUtils/pbarpStatesLS.hh"
#include "Examples/pbarpToOmegaPiLS/pbarpToOmegaPi0StatesLS.hh"
#include "ErrLogger/ErrLogger.hh"

AbsOmegaPiLhLS::AbsOmegaPiLhLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates) :
  _globalItCounter(0),
  _omegaPiEventListPtr(theEvtList),
  _omegaPi0StatesPtr(theStates)
{
  _evtDataVec=_omegaPiEventListPtr->getDataVecs();
  _evtMCVec=_omegaPiEventListPtr->getMcVecs();

  _allJPCLSlsStates=theStates->all_JPCLSls_States();
  _singlet_JPCLS_States=theStates->singlet_JPCLSls_States();
  _triplet0_JPCLS_States=theStates->triplet0_JPCLSls_States();
  _tripletp1_JPCLS_States=theStates->tripletp1_JPCLSls_States();
  _tripletm1_JPCLS_States=theStates->tripletm1_JPCLSls_States();
}

AbsOmegaPiLhLS::AbsOmegaPiLhLS(boost::shared_ptr<AbsOmegaPiLhLS> theAbsOmegaPiLhLSPtr):
  _omegaPiEventListPtr(theAbsOmegaPiLhLSPtr->getEventList()),
  _omegaPi0StatesPtr(theAbsOmegaPiLhLSPtr->omegaPi0States())
{
  _evtDataVec=_omegaPiEventListPtr->getDataVecs();
  _evtMCVec=_omegaPiEventListPtr->getMcVecs();

  _allJPCLSlsStates=_omegaPi0StatesPtr->all_JPCLSls_States();
  _singlet_JPCLS_States=_omegaPi0StatesPtr->singlet_JPCLSls_States();
  _triplet0_JPCLS_States= _omegaPi0StatesPtr->triplet0_JPCLSls_States();
  _tripletp1_JPCLS_States= _omegaPi0StatesPtr->tripletp1_JPCLSls_States();
  _tripletm1_JPCLS_States= _omegaPi0StatesPtr->tripletm1_JPCLSls_States();
}

AbsOmegaPiLhLS::~AbsOmegaPiLhLS()
{
}

double AbsOmegaPiLhLS::calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal){

  _globalItCounter++;
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<OmegaPiDataLS::OmPiEvtDataLS*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*>::iterator iterm;
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

void AbsOmegaPiLhLS::print(std::ostream& os) const{
  os << "AbsOmegaPiLhLS::print\n";
}
