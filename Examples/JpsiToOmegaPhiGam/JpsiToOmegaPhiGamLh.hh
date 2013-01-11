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
#include "Examples/JpsiToOmegaPhiGam/JpsiToOmegaPhiGamStates.hh"

#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;

class JpsiToOmegaPhiGamLh : public PsiProdBaseLhNew{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiToOmegaPhiGamLh(boost::shared_ptr<const EvtDataBaseListNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToOmegaPhiGamStates> theStates);
  JpsiToOmegaPhiGamLh(boost::shared_ptr<AbsLhNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiToOmegaPhiGamStates> theStates);

  /** Destructor */
  virtual ~JpsiToOmegaPhiGamLh();

  virtual AbsLhNew* clone_() const {
    return new JpsiToOmegaPhiGamLh(_evtListPtr, _hypVec, _jpsiToOmegaPhiGamStatesPtr);
  }
  
//   virtual double calcEvtIntensity( EvtDataNew* theData, fitParamsNew& theParamVal);
  
  //Getters:
//   virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  
  
protected:
  boost::shared_ptr<JpsiToOmegaPhiGamStates> _jpsiToOmegaPhiGamStatesPtr; 

private:

  void initializeHypothesis();

};



