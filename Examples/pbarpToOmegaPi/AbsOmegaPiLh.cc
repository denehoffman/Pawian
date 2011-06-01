#include <getopt.h>
#include <fstream>
#include <string>

#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiEventList.hh"
#include "PwaUtils/pbarpStates.hh"
#include "Examples/pbarpToOmegaPi/pbarpToOmegaPi0States.hh"
#include "ErrLogger/ErrLogger.hh"

AbsOmegaPiLh::AbsOmegaPiLh(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates) :
  _globalItCounter(0),
  _omegaPiEventListPtr(theEvtList),
  _omegaPi0StatesPtr(theStates)
{
  _evtDataVec=_omegaPiEventListPtr->getDataVecs();
  _evtMCVec=_omegaPiEventListPtr->getMcVecs();
}

AbsOmegaPiLh::AbsOmegaPiLh(boost::shared_ptr<AbsOmegaPiLh> theAbsOmegaPiLhPtr):
  _omegaPiEventListPtr(theAbsOmegaPiLhPtr->getEventList()),
  _omegaPi0StatesPtr(theAbsOmegaPiLhPtr->omegaPi0States())
{
  _evtDataVec=_omegaPiEventListPtr->getDataVecs();
  _evtMCVec=_omegaPiEventListPtr->getMcVecs();
}

AbsOmegaPiLh::~AbsOmegaPiLh()
{
}

double AbsOmegaPiLh::calcLogLh(const OmegaPiData::fitParamVal& theParamVal){

  _globalItCounter++;
 
  double logLH=0.;
  double logLH_data=0.;

  std::vector<OmegaPiData::OmPiEvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=log10(intensity);
  } 

  double LH_mc=0.;
  
  std::vector<OmegaPiData::OmPiEvtData*>::iterator iterm;
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

void AbsOmegaPiLh::print(std::ostream& os) const{
  os << "AbsOmegaPiLh::print\n";
}
