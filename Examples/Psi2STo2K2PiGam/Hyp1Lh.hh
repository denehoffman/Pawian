#ifndef _Hyp1Lh_H
#define _Hyp1Lh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <map>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Examples/Psi2STo2K2PiGam/AbsPsi2STo2K2PiGamLh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Hyp1Lh : public AbsPsi2STo2K2PiGamLh{

public:

  // create/copy/destroy:

  ///Constructor 
  Hyp1Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, const std::map<const std::string, bool>& hypMap);
  Hyp1Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Hyp1Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_(){
    Hyp1Lh* result = new Hyp1Lh(_Psi2STo2K2PiGamEvtListPtr, _hypMap);
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

  bool _K1_1270Hyp;
  bool _K1_1400Hyp;
  bool _K0_1430_K0_1430Hyp;
  bool _K2_1430_K2_1430Hyp;
  bool _K0_1430_K2_1430Hyp;
  bool _K0_1430_K892Hyp1;
  bool _K2_1430_K892Hyp1;
  bool _K1_1410_K1_1410Hyp;
  bool _K1_1410_K892Hyp1;
  bool _f1710_f1710Hyp1;
  std::map<const std::string, bool> _hypMap;

  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);
  virtual bool equalChic0DecParams();

  virtual void copyCurrentVals(Hyp1Lh* theLh);
  std::map<unsigned int, complex<double> > _currentResultHyp1;
  std::map<unsigned int, complex<double> > _currentResultlK2_1430_K2_1430Hyp1;
  std::map<unsigned int, complex<double> > _currentResultlK892K892Hyp1;
  std::map<unsigned int, complex<double> > _currentResultlK2_1430_K892Hyp1;

private:
  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;
  std::vector<unsigned int> _ampVecComp;
  std::vector<unsigned int> _ampVecK2_1430_K2_1430;
  std::vector<unsigned int> _massVecK2_1430_K2_1430;
  std::vector<unsigned int> _ampVecK892K892;
  std::vector<unsigned int> _massVecK892K892;
  std::vector<unsigned int> _ampVecK2_1430_K892;
  std::vector<unsigned int> _massVecK2_1430_K892;

  void setUp(const std::map<const std::string, bool>& hypMap);
  unsigned int _evtCounter;
  bool _equalParameter;
  bool _equalK2_1430_K2_1430Params;
  bool _equalK892K892Params;
  bool _equalK2_1430_K892Params;
  bool equalParams();
  
};

#endif
