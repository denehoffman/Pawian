#ifndef _Eta2DecAmp_H
#define _Eta2DecAmp_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"



class Eta2DecAmp : public XDecAmpBase{

public:

  // create/copy/destroy:

  ///Constructor 
  Eta2DecAmp(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);

  /** Destructor */
  virtual ~Eta2DecAmp();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;

protected:

private:

//   complex<double> EtaToAPiFlatteAmp(EvtData* theData, fitParamsNew& theParamVal);

};

#endif


