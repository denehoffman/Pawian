#ifndef _Hyp9Lh_H
#define _Hyp9Lh_H

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

#include "Examples/Psi2STo2K2PiGam/Hyp8Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Hyp9Lh : public Hyp8Lh{

public:

  // create/copy/destroy:

  ///Constructor 
  Hyp9Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, const std::map<const std::string, bool>& hypMap);
  Hyp9Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Hyp9Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() const{
    return new Hyp9Lh(_Psi2STo2K2PiGamEvtListPtr, _hypMap);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal);
  virtual int setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par);
  virtual unsigned int nFitParams();

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(param2K2PiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix);

protected:
  bool _Pi_2_1670Tof_2_1270PiHyp9;
  bool _Pi_2_1670ToK892KHyp9;
  bool _Pi1800Tof980PiHyp9;
  bool _doHyp9;

  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);

private:

  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;

  void setUp(const std::map<const std::string, bool>& hypMap);
};

#endif
