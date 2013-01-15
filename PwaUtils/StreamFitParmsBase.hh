#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsFitParamStreamer.hh"
#include "PwaUtils/FitParamsBase.hh"

class AbsLh;

class StreamFitParmsBase : public AbsFitParamStreamer {
  
public:
  StreamFitParmsBase(std::string&, boost::shared_ptr<AbsLh> theLhPtr);
  virtual ~StreamFitParmsBase();
  
  fitParams getFitParamVal() { return _paramVal;}
  fitParams getFitParamErr() { return _paramErr;}
 
  virtual void fillParamMap() {return;}

protected:
  virtual void fillParams();


  virtual void fillLamLamAmps(mapStrJPCLamLam& valMap, mapStrJPCLamLam& errMap, 
			      const std::string& suffix);
  virtual void fillLSAmps(mapStrJPCLS& valMap, mapStrJPCLS& errMap, 
			  const std::string& suffix);
  virtual void fillDoubles(mapStrDouble& valMap, mapStrDouble& errMap, 
			   const std::string& suffix);
private:
  fitParams _paramVal;
  fitParams _paramErr;

  void fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, 
		     std::string& suffix, int index);
};
