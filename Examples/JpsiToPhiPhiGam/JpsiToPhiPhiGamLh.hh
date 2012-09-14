#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "PwaUtils/PsiProdBaseLhNew.hh"
#include "Examples/JpsiToPhiPhiGam/JpsiToPhiPhiGamStates.hh"

#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;

class JpsiToPhiPhiGamLh : public PsiProdBaseLhNew{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiToPhiPhiGamLh(boost::shared_ptr<const EvtDataBaseListNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToPhiPhiGamStates> theStates);
  JpsiToPhiPhiGamLh(boost::shared_ptr<AbsLhNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToPhiPhiGamStates> theStates);

  /** Destructor */
  virtual ~JpsiToPhiPhiGamLh();

  virtual AbsLhNew* clone_() const {
    return new JpsiToPhiPhiGamLh(_evtListPtr, _hypVec, _jpsiToPhiPhiGamStatesPtr);
  }
  
//   virtual double calcEvtIntensity( EvtDataNew* theData, fitParamsNew& theParamVal);
  
  //Getters:
//   virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  
  
protected:
  boost::shared_ptr<JpsiToPhiPhiGamStates> _jpsiToPhiPhiGamStatesPtr; 


private:

  void initializeHypothesis();

};



