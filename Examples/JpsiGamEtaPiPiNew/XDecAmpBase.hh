#ifndef _XDecAmpBase_H
#define _XDecAmpBase_H

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"


class XDecAmpBase : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  XDecAmpBase(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates, Spin spinX);

  /** Destructor */
  virtual ~XDecAmpBase();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal)=0;
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr)=0;
  virtual void print(std::ostream& os) const;

protected:
  const std::string _a980piKey;
  bool _a980piHyp;
  const std::string _sigmaEtaKey;
  bool _sigmaEtaHyp;
  const double _massPi0;
  const double _massKplus;
  const double _massK0;
  const double _massEta;
  std::pair <const double, const double> _decPairKK;
  std::pair <const double, const double> _decPairPi0Eta;
  boost::shared_ptr<JpsiGamEtaPiPiStates> _theStatesPtr;
  Spin _spinX; 

  complex<double> XToAPiFlatteAmp(Spin lamX, EvtData* theData, fitParamsNew& theParamVal);
  complex<double> XToEtaFAmp(Spin lamX, Spin jf, EvtData* theData, 
			     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFMag, 
			     std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess >& XToEtaFPhi,
			     double fMass, double fWidth);
  void initialize();

private:




};

#endif


