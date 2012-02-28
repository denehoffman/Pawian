#ifndef _Hyp3Lh_H
#define _Hyp3Lh_H

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

#include "Examples/Psi2STo2K2PiGam/Hyp2Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Hyp3Lh : public Hyp2Lh{

public:

  // create/copy/destroy:

  ///Constructor 
  Hyp3Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, const std::map<const std::string, bool>& hypMap, bool chacheAmps=false);
  Hyp3Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap, bool chacheAmps=false);

  /** Destructor */
  virtual ~Hyp3Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_(){
    return new Hyp3Lh(_Psi2STo2K2PiGamEvtListPtr, _hypMap, _cacheAmps);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal);
  virtual int setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par);
  virtual unsigned int nFitParams();

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(param2K2PiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix);

protected:  
  bool _sigmaf980Hyp;
  bool _sigmaf1710Hyp;
  bool _sigmaf2200Hyp;
  bool _doHyp3;
  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);
  virtual bool equalChic0DecParams();

  virtual void copyCurrentVals(Hyp3Lh* theLh);
  std::map<unsigned int, complex<double> > _currentResultHyp3;
private:

  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;
  std::vector<unsigned int> _massVecRemain;
  bool _f980FlatteRemain;

  void setUp(const std::map<const std::string, bool>& hypMap);

  unsigned int _evtCounter;
  bool _equalParameter;

  bool equalParams();
};

#endif
