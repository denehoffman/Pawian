#ifndef Stream2K2PiGamFitParms_H
#define Stream2K2PiGamFitParms_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PwaUtils/AbsFitParamStreamer.hh"
#include "Examples/Psi2STo2K2PiGam/FitParams2K2PiGam.hh"


class Stream2K2PiGamFitParms : public AbsFitParamStreamer{
  
public:
  Stream2K2PiGamFitParms(std::string&);
  virtual ~Stream2K2PiGamFitParms();
  
  param2K2PiGam getFitParamVal() { return _paramVal;}
  param2K2PiGam getFitParamErr() { return _paramErr;}
  
protected:
  virtual void fillParamMap();

private:
  param2K2PiGam _paramVal;
  param2K2PiGam _paramErr;

  FitParams2K2PiGam _fitParams2K2PiGam;
};

#endif 
