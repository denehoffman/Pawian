#ifndef _Hyp5Lh_H
#define _Hyp5Lh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Examples/Psi2STo2K2PiGam/Hyp4Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Hyp5Lh : public Hyp4Lh{

public:

  // create/copy/destroy:

  ///Constructor 
  Hyp5Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, boost::shared_ptr<const Psi2STo2K2PiGamStates>, bool K1_1270Hyp, bool K0_1430_K0_1430Hyp, bool disableHyp3);
  Hyp5Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, bool K1_1270Hyp, bool K0_1430_K0_1430Hyp, bool disableHyp3);

  /** Destructor */
  virtual ~Hyp5Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() const {
    return new Hyp5Lh(_Psi2STo2K2PiGamEvtListPtr, _Psi2STo2K2PiGamStatesPtr, _K1_1270Hyp, _K0_1430_K0_1430Hyp, _disableHyp3);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal);
  virtual int setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const;
  virtual void dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const;

protected:
  virtual complex<double> chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

private:

};

#endif
