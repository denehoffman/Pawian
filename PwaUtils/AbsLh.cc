#include <getopt.h>
#include <fstream>
#include <string>

#include "PwaUtils/AbsLh.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

AbsLh::AbsLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  _evtListPtr(theEvtList)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
}

AbsLh::AbsLh(boost::shared_ptr<AbsLh> theAbsLhPtr):
  _evtListPtr(theAbsLhPtr->getEventList())
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
}

AbsLh::~AbsLh()
{
}

double AbsLh::calcLogLh(fitParams& theParamVal){
 
  double logLH=0.;
  double logLH_data=0.;
  double weightSum=0.;

  std::vector<EvtData*>::iterator iterd;
  for (iterd=_evtDataVec.begin(); iterd!=_evtDataVec.end(); ++iterd){
    double intensity=calcEvtIntensity((*iterd), theParamVal);
    if (intensity>0.) logLH_data+=((*iterd)->evtWeight)*log(intensity);
    weightSum+= (*iterd)->evtWeight;
  } 

  double LH_mc=0.;
  
  std::vector<EvtData*>::iterator iterm;
  for (iterm=_evtMCVec.begin(); iterm!=_evtMCVec.end(); ++iterm){
           double intensity=calcEvtIntensity((*iterm), theParamVal);
           LH_mc+=intensity;
         }

  double logLH_mc_Norm=0.;
  if (LH_mc>0.) logLH_mc_Norm=log(LH_mc/_evtMCVec.size());

  logLH=0.5*weightSum *(LH_mc/_evtMCVec.size()-1.)*(LH_mc/_evtMCVec.size()-1.)
    -logLH_data
    +weightSum*logLH_mc_Norm;

  Info << "current LH = " << logLH << endmsg;

 return logLH;

}

void AbsLh::setHyps( const std::map<const std::string, bool>& theMap, bool& theHyp, std::string& theKey){

  std::map<const std::string, bool>::const_iterator iter= theMap.find(theKey);
  
  if (iter !=theMap.end()){
    theHyp= iter->second;
    DebugMsg<< "hypothesis " << iter->first << "\t" << theHyp <<endmsg;
    _hypMap[iter->first]= iter->second;
  }
  else{
    Alert << theKey << " does not exist!!!" <<endmsg;
    exit(0);
  }
}

