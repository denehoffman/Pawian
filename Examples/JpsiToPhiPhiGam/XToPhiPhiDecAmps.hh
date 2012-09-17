#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamStates.hh"


class XToPhiPhiDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  XToPhiPhiDecAmps(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToPhiPhiGamStates> theStates, Spin spinX, int parity=1);

  /** Destructor */
  virtual ~XToPhiPhiDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;

protected:
  const std::string _phiPhiKey;
  const std::string _xBWKey;

  bool _massIndependent;

  boost::shared_ptr<JpsiToPhiPhiGamStates> _theStatesPtr;

  complex<double> XToPhiPhiAmp(Spin lamX, EvtDataNew* theData, fitParamsNew& theParamVal);
  complex<double> phiphiTo4KAmp( EvtDataNew* theData, Spin lambdaPhi1, Spin lambdaPhi2  );

  virtual void initialize();

private:




};



