#ifndef _Hyp7Lh_H
#define _Hyp7Lh_H

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

#include "Examples/Psi2STo2K2PiGam/Hyp6Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Hyp7Lh : public Hyp6Lh{

public:

  // create/copy/destroy:

  ///Constructor 
  Hyp7Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, boost::shared_ptr<const Psi2STo2K2PiGamStates>, const std::map<const std::string, bool>& hypMap);
  Hyp7Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Hyp7Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() const{
    return new Hyp7Lh(_Psi2STo2K2PiGamEvtListPtr, _Psi2STo2K2PiGamStatesPtr, _hypMap);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, Psi2STo2K2PiGamData::fitParamVal& startVal,  Psi2STo2K2PiGamData::fitParamVal& errVal);
  virtual int setFitParamVal(Psi2STo2K2PiGamData::fitParamVal& theParamVal, const std::vector<double>& par) const;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(Psi2STo2K2PiGamData::fitParamVal& theParamVal) const;
  virtual void dumpCurrentResult(std::ostream& os, Psi2STo2K2PiGamData::fitParamVal& theParamVal, std::string& suffix) const;

protected:

  bool _K1_1680Hyp;
  bool _K1_2300Hyp;

  virtual complex<double> chi0DecAmps(const Psi2STo2K2PiGamData::fitParamVal& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

private:

};

#endif
