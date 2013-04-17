//************************************************************************//
//									  //
//  Copyright 2013 Bertram Kopf (bertram@ep1.rub.de)			  //
//  	      	   Julian Pychy (julian@ep1.rub.de)			  //
//          	   - Ruhr-Universität Bochum 				  //
//									  //
//  This file is part of Pawian.					  //
//									  //
//  Pawian is free software: you can redistribute it and/or modify	  //
//  it under the terms of the GNU General Public License as published by  //
//  the Free Software Foundation, either version 3 of the License, or 	  //
//  (at your option) any later version.	 	      	  	   	  //
//									  //
//  Pawian is distributed in the hope that it will be useful,		  //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of	  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	  //
//  GNU General Public License for more details.	      		  //
//									  //
//  You should have received a copy of the GNU General Public License     //
//  along with Pawian.  If not, see <http://www.gnu.org/licenses/>.	  //
//									  //
//************************************************************************//

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
#include "PwaUtils/AbsEnv.hh"
#include "PwaUtils/ParserBase.hh"
#include "qft++/relativistic-quantum-mechanics/Utils.hh"
#include "ErrLogger/ErrLogger.hh"


AbsLh::AbsLh(boost::shared_ptr<AbsLh> theAbsLhPtr):
  AbsParamHandler()
  ,_absEnv(theAbsLhPtr->_absEnv)
  ,_evtDataVec(theAbsLhPtr->getDataVec())
  ,_evtMCVec(theAbsLhPtr->getMcVec())
  ,_usePhasespace(theAbsLhPtr->_absEnv->parser()->usePhaseSpaceHyp())
  ,_phasespaceKey("Phasespace")
  ,_calcCounter(0)
{
  unsigned int noOfThreads=_absEnv->parser()->noOfThreads();
  if(noOfThreads > boost::thread::hardware_concurrency()) noOfThreads=boost::thread::hardware_concurrency();
  _noOfThreads = noOfThreads;
}

AbsLh::AbsLh(AbsEnv* theEnv) :
  AbsParamHandler()
  ,_absEnv(theEnv)
  ,_usePhasespace(theEnv->parser()->usePhaseSpaceHyp())
  ,_phasespaceKey("Phasespace")
  ,_calcCounter(0)
{
  //  _noOfThreads = boost::thread::hardware_concurrency();
  _noOfThreads = theEnv->parser()->noOfThreads();
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

  LHData theLHData;
  theLHData.logLH_data = theLHData.weightSum = theLHData.LH_mc = 0.0;

  int eventStepData = _evtDataVec.size() / _noOfThreads;
  int eventStepMC = _evtMCVec.size() / _noOfThreads;

  std::vector<std::thread> theThreads;
  std::vector<LHData> threadDataVec;
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
     theLHData.logLH_data += (*it).logLH_data;
     theLHData.weightSum += (*it).weightSum;
     theLHData.LH_mc  += (*it).LH_mc;
  }

  return mergeLogLhData(theLHData, _evtMCVec.size());
}




void AbsLh::calcLogLhDataClient(fitParams& theParamVal,
				LHData& theLHData, std::vector<double>& eventLimits ){

  _calcCounter++;
  if (_cacheAmps && _calcCounter>1) checkRecalculation(theParamVal);
  updateFitParams(theParamVal);

  int numData = _evtDataVec.size();
  int numMC = _evtMCVec.size();

  int eventStepData = numData / _noOfThreads;
  int eventStepMC = numMC / _noOfThreads;

  std::vector<std::thread> theThreads;
  std::vector<LHData> threadDataVec;
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
     theLHData.logLH_data += (*it).logLH_data;
     theLHData.weightSum += (*it).weightSum;
     theLHData.LH_mc += (*it).LH_mc;
  }
}



double AbsLh::mergeLogLhData(LHData& theLHData, int nMCs){

  double logLH=0.;
  double logLH_mc_Norm=0.;  

  // if (theLHData.LH_mc>0.) logLH_mc_Norm=log(theLHData.LH_mc/_evtMCVec.size());
  // logLH=0.5*theLHData.weightSum *(theLHData.LH_mc/_evtMCVec.size()-1.)*(theLHData.LH_mc/_evtMCVec.size()-1.)
  //   -theLHData.logLH_data
  //   +theLHData.weightSum*logLH_mc_Norm;

  if (theLHData.LH_mc>0.) logLH_mc_Norm=log(theLHData.LH_mc/nMCs);
  logLH=0.5*theLHData.weightSum *(theLHData.LH_mc/nMCs-1.)*(theLHData.LH_mc/nMCs-1.)
    -theLHData.logLH_data
    +theLHData.weightSum*logLH_mc_Norm;  

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

  if(_usePhasespace){
    fitVal.otherParams[_phasespaceKey]=0.01;
    fitErr.otherParams[_phasespaceKey]=0.05;
  } 

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

void AbsLh::setDataVec(std::vector<EvtData*> theVec) {
  if(_evtDataVec.size()>0){
    Alert << "data vector already set!!!" << endmsg;
    exit(0); 
  }

  _evtDataVec=theVec;
}

void AbsLh::setMcVec(std::vector<EvtData*> theVec) {
  if(_evtMCVec.size()>0){
    Alert << "mc vector already set!!!" << endmsg;
    exit(0); 
  }

  _evtMCVec=theVec;
}



