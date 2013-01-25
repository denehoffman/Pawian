#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"
#include "Examples/JpsiToOmegaPhiGam/JpsiToOmegaPhiGamStates.hh"


class XToOmegaPhiDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  XToOmegaPhiDecAmps(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToOmegaPhiGamStates> theStates, Spin spinX, int parity=1);

  /** Destructor */
  virtual ~XToOmegaPhiDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParamsNew& theParamVal);

  void updateFitParams(fitParamsNew& theParamVal);
protected:
  const std::string _phiPhiKey;
  const std::string _xBWKey;
  const std::string _xFlatteKey;

  bool _massIndependent;
  bool _bwMassFit;
  bool _flatteMassFit;

  const double _phiMass;
  const double _omegaMass;
  pair<const double, const double> _phiPhiPair;
  pair<const double, const double> _omegaPhiPair;

  boost::shared_ptr<JpsiToOmegaPhiGamStates> _theStatesPtr;

  complex<double> XToPhiOmegaAmp(Spin lamX, EvtDataNew* theData);
  complex<double> phiOmegaTo2K3piAmp( EvtDataNew* theData, Spin lambdaOmega, Spin lambdaPhi  );

  virtual void initialize();

private:




};



