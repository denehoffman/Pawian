#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/AbsLhNew.hh"
#include "Examples/D0ToKsPipPim/D0ToKsPipPimStates.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"
#include "PwaDynamics/FVector.hh"
#include "PwaDynamics/FVectorPiPiS.hh"

class AbsXdecAmp;

class D0ToKsPipPimLh : public AbsLhNew {

public:
  D0ToKsPipPimLh(boost::shared_ptr<const EvtDataBaseListNew>, 
		   const std::vector<std::string>& hypVec, 
		   boost::shared_ptr<D0ToKsPipPimStates> theStates);
  D0ToKsPipPimLh(boost::shared_ptr<AbsLhNew>, const std::vector<std::string>& hypVec, 
		   boost::shared_ptr<D0ToKsPipPimStates> theStates);
  virtual ~D0ToKsPipPimLh();

  virtual AbsLhNew* clone_() const{
    return new D0ToKsPipPimLh(_evtListPtr, _hypVec, _theStatesPtr);
  }
  virtual double calcEvtIntensity( EvtDataNew* theData, fitParamsNew& theParamVal);
  
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;

protected:
  const std::vector<std::string> _hypVec;
  boost::shared_ptr<D0ToKsPipPimStates> _theStatesPtr;
 
  bool _usePhasespace;
  const std::string _phasespaceKey;

  bool _pipiSHyp;
  const std::string _pipiSKey;

  bool _KpiSHyp;
  const std::string _KpiSKey;

  std::map< std::string,std::vector<std::string> > _hypMaps;

private:

  void getDefaultLSParams(const std::vector<std::string>& hyps, 
			  std::vector< boost::shared_ptr<const JPCLS> > lsAmps, 
			  fitParamsNew& fitVal, fitParamsNew& fitErr);
  void initializeHypothesis();


};
