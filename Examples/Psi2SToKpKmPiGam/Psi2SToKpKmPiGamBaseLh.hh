#ifndef _Psi2SToKpKmPiGamBaseLh_H
#define _Psi2SToKpKmPiGamBaseLh_H

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

#include "Examples/Psi2SToKpKmPiGam/AbsPsi2SToKpKmPiGamLh.hh"
#include "Examples/Psi2SToKpKmPiGam/Psi2SToKpKmPiGamData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class Psi2SToKpKmPiGamBaseLh : public AbsPsi2SToKpKmPiGamLh{

public:

  // create/copy/destroy:

  ///Constructor 
  Psi2SToKpKmPiGamBaseLh(boost::shared_ptr<const Psi2SToKpKmPiGamEventList>, const std::map<const std::string, bool>& hypMap);
  Psi2SToKpKmPiGamBaseLh(boost::shared_ptr<AbsPsi2SToKpKmPiGamLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~Psi2SToKpKmPiGamBaseLh();

  virtual AbsPsi2SToKpKmPiGamLh* clone_() const {
    return new Psi2SToKpKmPiGamBaseLh(_Psi2SToKpKmPiGamEvtListPtr, _hypMap);
  }


  // Getters:
  virtual void setMnUsrParams(MnUserParameters& upar, paramKpKmPiGam& startVal,  paramKpKmPiGam& errVal);
  virtual int setFitParamVal(paramKpKmPiGam& theParamVal, const std::vector<double>& par);

  virtual unsigned int nFitParams();
  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(paramKpKmPiGam& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, paramKpKmPiGam& theParamVal, std::string& suffix);
protected:
  bool _K0_1430Hyp;
  bool _K1_1410Hyp;
  bool _K2_1430Hyp;
  bool _KKPi_Hyp;
  std::map<const std::string, bool> _hypMap;
  virtual  complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData);

  virtual  complex<double> calcDecAmp(complex<double>& inAmp, Spin lamChi, const paramKpKmPiGam& theParamVal, Psi2SToKpKmPiGamData::Psi2SToKpKmPiGamEvtData* theData);
private:
  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;


};

#endif
