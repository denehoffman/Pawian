#ifndef _JpsiGamEtaPiPiProdLhNew_H
#define _JpsiGamEtaPiPiProdLhNew_H

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
#include "Examples/JpsiGamEtaPiPi/JpsiGamEtaPiPiStates.hh"

#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;

class JpsiGamEtaPiPiProdLhNew : public PsiProdBaseLhNew{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamEtaPiPiProdLhNew(boost::shared_ptr<const EvtDataBaseListNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);
  JpsiGamEtaPiPiProdLhNew(boost::shared_ptr<AbsLhNew>, const std::vector<std::string>& hypVec, boost::shared_ptr<JpsiGamEtaPiPiStates> theStates);

  /** Destructor */
  virtual ~JpsiGamEtaPiPiProdLhNew();

  virtual AbsLhNew* clone_() const {
    return new JpsiGamEtaPiPiProdLhNew(_evtListPtr, _hypVec, _jpsiGamEtaPiPiStatesPtr);
  }
  
//   virtual double calcEvtIntensity( EvtDataNew* theData, fitParamsNew& theParamVal);
  
  //Getters:
//   virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  
  
protected:
  boost::shared_ptr<JpsiGamEtaPiPiStates> _jpsiGamEtaPiPiStatesPtr; 


private:

  void initializeHypothesis();


  

};

#endif


