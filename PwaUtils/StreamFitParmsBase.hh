#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsFitParamStreamer.hh"
#include "PwaUtils/FitParamsBase.hh"

class StreamFitParmsBase : public AbsFitParamStreamer {
  
public:
  StreamFitParmsBase(std::string&, boost::shared_ptr<FitParamsBase> fitParamsBasePtr);
  virtual ~StreamFitParmsBase();
  
  fitParams getFitParamVal() { return _paramVal; }
  fitParams getFitParamErr() { return _paramErr; }
   virtual void fillParamMap() { return; }

protected:
  virtual void fillParams();

  virtual void fillAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, 
			std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, 
			double, pawian::Collection::SharedPtrLess >& valMap , 
			std::map< boost::shared_ptr<const JPCLS>, double,
			pawian::Collection::SharedPtrLess >& errMap);

  virtual void fillLamLamAmps(std::vector< boost::shared_ptr<const JPClamlam> >& theJPCLamLams,
			      std::string& suffix, std::map< boost::shared_ptr<const JPClamlam>,
			      double, pawian::Collection::SharedPtrLess >& valMap, 
			      std::map< boost::shared_ptr<const JPClamlam>, double, 
			      pawian::Collection::SharedPtrLess >& errMap);
  
private:
  fitParams _paramVal;
  fitParams _paramErr;
  boost::shared_ptr<FitParamsBase> _fitParamsBasePtr;

  void fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap,
		     std::string& suffix, int index);
};
