#ifndef _Hyp4Lh_H
#define _Hyp4Lh_H

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

#include "Examples/Psi2STo2K2PiGam/Hyp3Lh.hh"
#include "Examples/Psi2STo2K2PiGam/Psi2STo2K2PiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Hyp4Lh : public Hyp3Lh{

public:

  // create/copy/destroy:

  ///Constructor 
  Hyp4Lh(boost::shared_ptr<const Psi2STo2K2PiGamEvtList>, const std::map<const std::string, bool>& hypMap);
  Hyp4Lh(boost::shared_ptr<AbsPsi2STo2K2PiGamLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Hyp4Lh();

  virtual AbsPsi2STo2K2PiGamLh* clone_(){
    return new Hyp4Lh(_Psi2STo2K2PiGamEvtListPtr, _hypMap);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, param2K2PiGam& startVal,  param2K2PiGam& errVal);
  virtual int setFitParamVal(param2K2PiGam& theParamVal, const std::vector<double>& par);
  virtual unsigned int nFitParams();

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(param2K2PiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, param2K2PiGam& theParamVal, std::string& suffix);

protected:
  bool _f980f1370Hyp4;
  bool _f980f1500Hyp4;
  bool _sigmaf1370Hyp4;
  bool _f1710f1370Hyp4;
  bool _f980f_2_1270Hyp4;
  bool _f980f_2_1430Hyp4;
  bool _f980f_2_1525Hyp4;
  bool _f980f_2_1950Hyp4;
  bool _f980f_2_2010Hyp4;
  bool _f1500f_2_1525Hyp4;
  bool _f1710f_2_1430Hyp4;
  bool _f1710f_2_1950Hyp4;
  bool _doHyp4;
  virtual complex<double> chi0DecAmps(const param2K2PiGam& theParamVal, Psi2STo2K2PiGamData::Psi2STo2K2PiGamEvtData* theData);
  virtual bool equalChic0DecParams();

  virtual void copyCurrentVals(Hyp4Lh* theLh);
  std::map<unsigned int, complex<double> > _currentResultHyp4;

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
