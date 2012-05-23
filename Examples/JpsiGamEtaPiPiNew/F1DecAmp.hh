#ifndef _F1DecAmp_H
#define _F1DecAmp_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "Examples/JpsiGamEtaPiPiNew/XDecAmpBase.hh"



class F1DecAmp : public XDecAmpBase{

public:

  // create/copy/destroy:

  ///Constructor 
  F1DecAmp(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);

  /** Destructor */
  virtual ~F1DecAmp();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;

protected:

private:


};

#endif


