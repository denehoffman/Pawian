#ifndef StreamFitParmsBaseNew_H
#define StreamFitParmsBaseNew_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PwaUtils/AbsFitParamStreamer.hh"
#include "PwaUtils/FitParamsBaseNew.hh"

class AbsLhNew;

class StreamFitParmsBaseNew : public AbsFitParamStreamer{
  
public:
  StreamFitParmsBaseNew(std::string&, boost::shared_ptr<AbsLhNew> theLhPtr);
  virtual ~StreamFitParmsBaseNew();
  
  fitParamsNew getFitParamVal() { return _paramVal;}
  fitParamsNew getFitParamErr() { return _paramErr;}
 
  virtual void fillParamMap() {return;}

protected:
  virtual void fillParams();


  virtual void fillLamLamAmps(mapStrJPCLamLam& valMap, mapStrJPCLamLam& errMap, const std::string& suffix);
  virtual void fillLSAmps(mapStrJPCLS& valMap, mapStrJPCLS& errMap, const std::string& suffix);
  virtual void fillDoubles(mapStrDouble& valMap, mapStrDouble& errMap, const std::string& suffix);
private:
  fitParamsNew _paramVal;
  fitParamsNew _paramErr;


  void fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index);
};


#endif 
