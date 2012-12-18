#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimStates.hh"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/KMatrixKPiSFocus.hh"
#include "PwaDynamics/PVectorKPiSFocus.hh"

class D0ToKPiSPiDec : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
 D0ToKPiSPiDec(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates);

  /** Destructor */
  virtual ~D0ToKPiSPiDec();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParamsNew& theParamVal);

  virtual void updateFitParams(fitParamsNew& theParamVal);

protected:
  const std::string _iso1Key;
  const std::string _iso3Key;
  bool _iso1Hyp;
  bool _iso3Hyp;

  boost::shared_ptr<D0ToKsPipPimStates> _theStatesPtr;

  virtual void initialize();

private:
  boost::shared_ptr<FVector> _fVec;  
  boost::shared_ptr<KMatrixKPiSFocus> _kMatr;
  boost::shared_ptr<PVectorKPiSFocus> _pVec;

  std::map<std::string, double> _currentbFactorMap;
  std::map<std::string, double> _currentaProdMap; 
  std::map<std::string, double> _currentbProdMap; 
  std::map<std::string, double> _currentcProdMap;
  std::map<std::string, double> _currentphaseProdMap;  
};



