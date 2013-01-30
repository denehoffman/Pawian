// AbsLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"
#include "PwaUtils/AbsXdecAmp.hh"

class AbsLh {

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
  void cacheAmplitudes();
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr)=0;
  virtual void print(std::ostream& os) const=0;
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  boost::shared_ptr<const EvtDataBaseList> _evtListPtr;
  std::vector<EvtData*> _evtDataVec;
  std::vector<EvtData*> _evtMCVec;
  std::map<const std::string, bool> _hypMap;
  bool _cacheAmps;
  unsigned int _calcCounter;

  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _allDecAmpMap; 

  virtual void setHyps( const std::map<const std::string, bool>& theMap, 
			bool& theHyp, std::string& theKey);
  
  virtual void checkParamVariation(fitParams& theParamVal);
  // virtual void cacheTheAmps()=0;
  virtual void cacheTheAmps();
};
