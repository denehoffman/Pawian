#ifndef _EtaDecAmp_H
#define _EtaDecAmp_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"



class EtaDecAmp : public XDecAmpBase{

public:

  // create/copy/destroy:

  ///Constructor 
  EtaDecAmp(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);

  /** Destructor */
  virtual ~EtaDecAmp();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;

protected:

private:

//   complex<double> EtaToAPiFlatteAmp(EvtData* theData, fitParamsNew& theParamVal);

};

#endif


