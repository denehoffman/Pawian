#ifndef StreamChic1ToKpKmPiGamFitParms_H
#define StreamChic1ToKpKmPiGamFitParms_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PwaUtils/AbsFitParamStreamer.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsChic1ToKpKmPiGam.hh"


class StreamChic1ToKpKmPiGamFitParms : public AbsFitParamStreamer{
  
public:
  StreamChic1ToKpKmPiGamFitParms(std::string&);
  virtual ~StreamChic1ToKpKmPiGamFitParms();
  
  fitParams getFitParamVal() { return _paramVal;}
  fitParams getFitParamErr() { return _paramErr;}
 
  virtual void fillParamMap() {return;}

protected:
  virtual void fillParams();
  virtual void fillAmps(std::vector< boost::shared_ptr<const JPCLS> >& theJPCLSs, std::string& suffix, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& valMap , std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& errMap);
private:
  fitParams _paramVal;
  fitParams _paramErr;

  FitParamsChic1ToKpKmPiGam _fitParamsKpKmPiGam;

  void fillAmps(std::string& suffix, std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& theAmpMapVal, std::map<int, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > >& theAmpMapErr);
  void fillParameter(std::map<int, double>& theValMap, std::map<int, double>& theErrMap, std::string& suffix, int index);
};


#endif 
