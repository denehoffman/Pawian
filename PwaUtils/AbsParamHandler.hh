// AbsParamHandler class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseList.hh"
#include "PwaUtils/FitParamsBase.hh"

class AbsParamHandler {

public:
  AbsParamHandler();
  virtual ~AbsParamHandler();
  
  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr)=0;
  virtual bool checkRecalculation(fitParams& theParamVal)=0;
  virtual void cacheAmplitudes()=0;
  virtual void updateFitParams(fitParams& theParamVal)=0;

protected:
  bool _cacheAmps;
  bool _recalculate;

private:

};
