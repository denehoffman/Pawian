#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/EvtDataBaseListNew.hh"
#include "PwaUtils/FitParamsBaseNew.hh"


class AbsLhNew {

public:
  AbsLhNew(boost::shared_ptr<const EvtDataBaseListNew>);
  AbsLhNew(boost::shared_ptr<AbsLhNew>);
  virtual ~AbsLhNew();
  virtual AbsLhNew* clone_() const = 0;

  virtual double calcLogLh(fitParamsNew& theParamVal);
  virtual double calcEvtIntensity(EvtDataNew* theData, fitParamsNew& theParamVal)=0;
  virtual boost::shared_ptr<const EvtDataBaseListNew> getEventList() const {
    return _evtListPtr;
  }
  void cacheAmplitudes();
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const=0;

protected:
  boost::shared_ptr<const EvtDataBaseListNew> _evtListPtr;
  std::vector<EvtDataNew*> _evtDataVec;
  std::vector<EvtDataNew*> _evtMCVec;
  std::map<const std::string, bool> _hypMap;
  bool _cacheAmps;

  virtual void setHyps( const std::map<const std::string, bool>& theMap, 
			bool& theHyp, std::string& theKey);
  
  virtual void checkParamVariation(fitParamsNew& theParamVal);
  virtual void cacheTheAmps()=0;
};
