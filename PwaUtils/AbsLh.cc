// AbsLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <iomanip>

#include "PwaUtils/AbsLh.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

AbsLh::AbsLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  AbsParamHandler()
  ,_evtListPtr(theEvtList)
  ,_calcCounter(0)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
}

AbsLh::AbsLh(boost::shared_ptr<AbsLh> theAbsLhPtr):
  AbsParamHandler()
  ,_evtListPtr(theAbsLhPtr->getEventList())
  ,_calcCounter(0)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
}

AbsLh::~AbsLh()
{
}

double AbsLh::calcLogLh(fitParams& theParamVal){
  _calcCounter++;
  if (_cacheAmps && _calcCounter>1) checkRecalculation(theParamVal); 
  updateFitParams(theParamVal);
  
  double logLH=0.;
  double logLH_data=0.;
  double weightSum=0.;
  double LH_mc=0.;


#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i=0; i<_evtDataVec.size(); ++i){
    EvtData* currentEvtData=_evtDataVec[i];
    double intensity=calcEvtIntensity(currentEvtData, theParamVal);

#ifdef _OPENMP
#pragma omp critical
      {
#endif
    if (intensity>0.) logLH_data+=(currentEvtData->evtWeight)*log(intensity);
    weightSum+= currentEvtData->evtWeight;

#ifdef _OPENMP
       }
#endif
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i=0; i<_evtMCVec.size(); ++i){
    EvtData* currentEvtData=_evtMCVec[i];
    double intensity=calcEvtIntensity(currentEvtData, theParamVal);

#ifdef _OPENMP
#pragma omp critical
      {
#endif
        LH_mc+=intensity;
#ifdef _OPENMP
       }
#endif
  }


  double logLH_mc_Norm=0.;  
  if (LH_mc>0.) logLH_mc_Norm=log(LH_mc/_evtMCVec.size());
  logLH=0.5*weightSum *(LH_mc/_evtMCVec.size()-1.)*(LH_mc/_evtMCVec.size()-1.)
    -logLH_data
    +weightSum*logLH_mc_Norm;
  
  Info << "current LH = " << std::setprecision(10) << logLH << endmsg;
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

void AbsLh::getDefaultParams(fitParams& fitVal, fitParams& fitErr){ 

  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->getDefaultParams(fitVal, fitErr);
  }
}

void AbsLh::cacheAmplitudes(){
  _cacheAmps=true;
  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->cacheAmplitudes();
  }
}

void AbsLh::updateFitParams(fitParams& theParamVal){

  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator itDecs;
  for(itDecs=_decAmps.begin(); itDecs!=_decAmps.end(); ++itDecs){
    (*itDecs)->updateFitParams(theParamVal);
  }
}

bool AbsLh::checkRecalculation(fitParams& theParamVal){
  bool result=true;
  std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    if(!(*it)->checkRecalculation(theParamVal)) result=false;
  }
  return result;
}


