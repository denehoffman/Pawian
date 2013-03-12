// AbsLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#include <getopt.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <thread>

#include <boost/thread.hpp>

#include "PwaUtils/AbsLh.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"


AbsLh::AbsLh(boost::shared_ptr<const EvtDataBaseList> theEvtList) :
  AbsParamHandler()
  ,_evtListPtr(theEvtList)
  ,_calcCounter(0)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
  _noOfThreads = boost::thread::hardware_concurrency();
}

AbsLh::AbsLh(boost::shared_ptr<AbsLh> theAbsLhPtr):
  AbsParamHandler()
  ,_evtListPtr(theAbsLhPtr->getEventList())
  ,_calcCounter(0)
{
  _evtDataVec=_evtListPtr->getDataVecs();
  _evtMCVec=_evtListPtr->getMcVecs();
  _noOfThreads = boost::thread::hardware_concurrency();
}

AbsLh::~AbsLh()
{
}



void AbsLh::ThreadfuncData(unsigned int minEvent, unsigned int maxEvent,
			   double& logLH_data, double& weightSum, fitParams& theParamVal){

   logLH_data=0.;
   weightSum=0.;

   for (unsigned int i=minEvent; i<=maxEvent; ++i){
      EvtData* currentEvtData=_evtDataVec[i];
      double intensity=calcEvtIntensity(currentEvtData, theParamVal);
      logLH_data+=(currentEvtData->evtWeight)*log(intensity);
      weightSum+= currentEvtData->evtWeight;
   }
}



void AbsLh::ThreadfuncMc(unsigned int minEvent, unsigned int maxEvent,
			 double& lh_mc, fitParams& theParamVal ){

   lh_mc=0.;

   for (unsigned int i=minEvent; i<=maxEvent; ++i){
      EvtData* currentEvtData=_evtMCVec[i];
      double intensity=calcEvtIntensity(currentEvtData, theParamVal);
      lh_mc+=intensity;
   }
}



double AbsLh::calcLogLh(fitParams& theParamVal){

  _calcCounter++;
  if (_cacheAmps && _calcCounter>1) checkRecalculation(theParamVal); 
  updateFitParams(theParamVal);

  double logLH=0.;
  double logLH_data=0.;
  double weightSum=0.;
  double LH_mc=0.;

  int eventStepData = _evtDataVec.size() / _noOfThreads;
  int eventStepMC = _evtMCVec.size() / _noOfThreads;

  std::vector<std::thread> theThreads;
  std::vector<ThreadData> threadDataVec;
  threadDataVec.resize(_noOfThreads);


  for(int i = 0; i<_noOfThreads;i++){

     int eventMin = i*eventStepData;
     int eventMax = (i==_noOfThreads-1) ? (_evtDataVec.size() - 1) : (i+1)*eventStepData - 1;

     theThreads.push_back(std::thread(&AbsLh::ThreadfuncData, this, eventMin, eventMax,
				      std::ref(threadDataVec.at(i).logLH_data), 
				      std::ref(threadDataVec.at(i).weightSum), theParamVal));
  }
  for(auto it = theThreads.begin(); it != theThreads.end(); ++it){
     (*it).join();
  }

  theThreads.clear();

  for(int i = 0; i<_noOfThreads;i++){

     int eventMin = i*eventStepMC;
     int eventMax = (i==_noOfThreads-1) ? (_evtMCVec.size() - 1) : (i+1)*eventStepMC - 1;

     theThreads.push_back(std::thread(&AbsLh::ThreadfuncMc, this, eventMin, eventMax,
				      std::ref(threadDataVec.at(i).LH_mc), theParamVal));
  }
  for(auto it = theThreads.begin(); it != theThreads.end(); ++it){
     (*it).join();
  }



  for(auto it = threadDataVec.begin(); it!= threadDataVec.end(); ++it){
     logLH_data += (*it).logLH_data;
     weightSum += (*it).weightSum;
     LH_mc += (*it).LH_mc;
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
std::vector< boost::shared_ptr<AbsXdecAmp> >::iterator it;
  for (it=_decAmps.begin(); it!=_decAmps.end(); ++it){
    (*it)->updateFitParams(theParamVal);
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


