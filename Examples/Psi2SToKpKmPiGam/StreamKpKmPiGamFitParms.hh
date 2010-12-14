#ifndef StreamKpKmPiGamFitParms_H
#define StreamKpKmPiGamFitParms_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PwaUtils/AbsFitParamStreamer.hh"
#include "Examples/Psi2SToKpKmPiGam/FitParamsKpKmPiGam.hh"


class StreamKpKmPiGamFitParms : public AbsFitParamStreamer{
  
public:
  StreamKpKmPiGamFitParms(std::string&);
  virtual ~StreamKpKmPiGamFitParms();
  
  paramKpKmPiGam getFitParamVal() { return _paramVal;}
  paramKpKmPiGam getFitParamErr() { return _paramErr;}
  
protected:
  virtual void fillParamMap();

private:
  paramKpKmPiGam _paramVal;
  paramKpKmPiGam _paramErr;

  FitParamsKpKmPiGam _fitParamsKpKmPiGam;
};

#endif 
