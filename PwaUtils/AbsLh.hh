// AbsLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsParamHandler.hh"
#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsXdecAmp.hh"


struct ThreadData{
   double weightSum;
   double logLH_data;
   double LH_mc;
};


class AbsLh : public AbsParamHandler{

public:
  AbsLh(boost::shared_ptr<const EvtDataBaseList>);
  AbsLh(boost::shared_ptr<AbsLh>);
  virtual ~AbsLh();
  virtual AbsLh* clone_() const = 0;

  virtual double calcLogLh(fitParams& theParamVal);
  virtual double calcEvtIntensity(EvtData* theData, fitParams& theParamVal)=0;
  virtual boost::shared_ptr<const EvtDataBaseList> getEventList() const {
    return _evtListPtr;
  }

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);  
  virtual void cacheAmplitudes();
  virtual void updateFitParams(fitParams& theParamVal);

  virtual void print(std::ostream& os) const=0;

protected:
  boost::shared_ptr<const EvtDataBaseList> _evtListPtr;
  std::vector<EvtData*> _evtDataVec;
  std::vector<EvtData*> _evtMCVec;
  std::vector< boost::shared_ptr<AbsXdecAmp> > _decAmps;

  std::map<const std::string, bool> _hypMap;
  unsigned int _calcCounter;
  unsigned short _noOfThreads;

  virtual void setHyps( const std::map<const std::string, bool>& theMap, 
			bool& theHyp, std::string& theKey);

  virtual void ThreadfuncData(unsigned int minEvent, unsigned int maxEvent,
			      double& logLH_data, double& weightSum, fitParams& theParamVal);
  virtual void ThreadfuncMc(unsigned int minEvent, unsigned int maxEvent,
			    double& lh_mc, fitParams& theParamVal);
};
