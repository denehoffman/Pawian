#ifndef Stream2K2PiGamFitParms_H
#define Stream2K2PiGamFitParms_H

#include <iostream>
#include <vector>
#include <map>
//#include <pair>
#include <string>
#include <boost/shared_ptr.hpp>
#include "PwaUtils/AbsFitParamStreamer.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamStates.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"


class Stream2K2PiGamFitParms : public AbsFitParamStreamer{
  
public:
  Stream2K2PiGamFitParms(std::string& , boost::shared_ptr<const Psi2STo2K2PiGamStates>);
  virtual ~Stream2K2PiGamFitParms();
  
  Psi2STo2K2PiGamData::fitParamVal getFitParamVal() { return _paramVal;}
  Psi2STo2K2PiGamData::fitParamVal getFitParamErr() { return _paramErr;}
  
protected:
  virtual void fillParamMap();

private:
  Psi2STo2K2PiGamData::fitParamVal _paramVal;
  Psi2STo2K2PiGamData::fitParamVal _paramErr;
  boost::shared_ptr<const Psi2STo2K2PiGamStates> _states;

};

#endif 
