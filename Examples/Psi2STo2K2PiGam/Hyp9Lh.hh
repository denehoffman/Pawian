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
  Hyp9Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps=false);
  Hyp9Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap, boost::shared_ptr<Psi2STo2K2PiGamStates> theStatesPtr, bool cacheAmps=false);

  /** Destructor */
  virtual ~Hyp9Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_() {
    Hyp9Lh* result = new Hyp9Lh(_Psi2STo2K2PiGamEvtListPtr, _hypMap, _fitParams2K2PiGam.states(), _cacheAmps);
    copyCurrentVals(result);
    return result;
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
  bool _Pi1800Tof1370PiHyp9;
  bool _Pi1800ToKappaKHyp9;
  bool _Pi1800Pi0ToK892KHyp9;
  bool _Pi3000Pi0ToK892KHyp9;
  bool _Pi3000Pi0ToK_0_1950KHyp9;
  bool _Pi_2_2285Tof1710PiHyp9;
  bool _Pi_2_2285ToK892KHyp9;
  bool _Pi_2_2285ToK_0_1430KHyp9;
  bool _Pi_2_2285ToK_2_1430KHyp9;
  bool _f980f_2_2300Hyp9;
  bool _f_2_2300sigmaHyp9;
  bool _K_2_1770ToK_2_1430PiHyp9;
  bool _doHyp9;

  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);
  virtual bool equalChic0DecParams();

  virtual void copyCurrentVals(Hyp9Lh* theLh);
  std::map<unsigned int, complex<double> > _currentResultHyp9;
  std::map<unsigned int, complex<double> > _currentResultPi_2_2285ToK_2_1430Hyp9;
  std::map<unsigned int, complex<double> > _currentResultPi_2_2285ToK892KHyp9;

private:

  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;
  std::vector<unsigned int> _massVecRemain;
  bool _f980FlatteRemain;

  std::vector<unsigned int> _ampVecPi_2_2285ToK_2_1430;
  std::vector<unsigned int> _massVecPi_2_2285ToK_2_1430;
  std::vector<unsigned int> _ampVecPi_2_2285ToK892K;
  std::vector<unsigned int> _massVecPi_2_2285ToK892K;

  void setUp(const std::map<const std::string, bool>& hypMap);

  unsigned int _evtCounter;
  bool _equalParameter;
  bool _equalPi_2_2285ToK_2_1430Params;
  bool _equalPi_2_2285ToK892KParams;

  bool equalParams();
};

#endif
