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
  Hyp7Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps=false);
  Hyp7Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr,bool cacheAmps=false);

  /** Destructor */
  virtual ~Hyp7Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() {
    return new Hyp7Lh(_Psi2STo2K2PiGamEvtListPtr, _hypMap, _fitParams2K2PiGam.states(), _cacheAmps);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal);
  virtual int setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par);
  virtual unsigned int nFitParams();

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(param2K2PiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix);

protected:
  bool _KappaHyp;
  bool _K1_1680Hyp;
  bool _K1_1680K1_1680Hyp7;
  bool _K1_1680K0_1430Hyp7;
  bool _K1_2300Hyp;

  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);
  virtual bool equalChic0DecParams();

  virtual void copyCurrentVals(Hyp7Lh* theLh);
  std::map<unsigned int, complex<double> > _currentResultHyp7;
  std::map<unsigned int, complex<double> > _currentResultK1_1680Hyp7;
  std::map<unsigned int, complex<double> > _currentResultK1_2300Hyp7;

private:

  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;
  std::vector<unsigned int> _massVecRemain;
  std::vector<unsigned int> _ampVecK1_1680;
  std::vector<unsigned int> _massVecK1_1680;
  std::vector<unsigned int> _ampVecK1_2300;
  std::vector<unsigned int> _massVecK1_2300;

  void setUp(const std::map<const std::string, bool>& hypMap);

  unsigned int _evtCounter;
  bool _equalParameter;
  bool _equalK1_1680Params;
  bool _equalK1_2300Params;

  bool equalParams();
};

#endif
