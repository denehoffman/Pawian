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
#include "PwaDynamics/FVectorPiPiS.hh"

class D0ToPiPiSKDec : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  D0ToPiPiSKDec(const std::string& name, const std::vector<std::string>& hypVec, boost::shared_ptr<D0ToKsPipPimStates> theStates);

  /** Destructor */
  virtual ~D0ToPiPiSKDec();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData);
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParamsNew& theParamVal);

  virtual void updateFitParams(fitParamsNew& theParamVal);

protected:
  const std::string _ASKey;
  bool _piPiSASHyp;
  boost::shared_ptr<D0ToKsPipPimStates> _theStatesPtr;

  virtual void initialize();

private:
  boost::shared_ptr<FVectorPiPiS> _pipiSFVec;  

  std::map<std::string, double> _currentgFactorMap;
  std::map<std::string, double> _currentbFactorMap;
  std::map<std::string, double> _currentfProdFactorMap; 
  double _currentS0Val;

  std::map<std::string, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > > _currentParamMagMap;
  std::map<std::string, std::map< boost::shared_ptr<const JPCLS>, double, pawian::Collection::SharedPtrLess > > _currentParamPhiMap;
};



