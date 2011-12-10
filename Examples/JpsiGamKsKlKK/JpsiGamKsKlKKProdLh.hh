#ifndef _JpsiGamKsKlKKProdLh_H
#define _JpsiGamKsKlKKProdLh_H

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

#include "Examples/JpsiGamKsKlKK/AbsJpsiGamKsKlKKLh.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"



class JpsiGamKsKlKKProdLh : public AbsJpsiGamKsKlKKLh{

public:

  // create/copy/destroy:

  ///Constructor 
  JpsiGamKsKlKKProdLh(boost::shared_ptr<const JpsiGamKsKlKKEventList>);
  JpsiGamKsKlKKProdLh(boost::shared_ptr<AbsJpsiGamKsKlKKLh>);

  /** Destructor */
  virtual ~JpsiGamKsKlKKProdLh();

  virtual AbsJpsiGamKsKlKKLh* clone_() const {
    return new JpsiGamKsKlKKProdLh(_JpsiGamKsKlKKEvtListPtr);
  }


  // Getters:
  virtual double calcEvtIntensity(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, const paramGamKsKlKK& theParamVal);
  
  virtual void setMnUsrParams(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal);
  virtual int setFitParamVal(paramGamKsKlKK& theParamVal, const std::vector<double>& par);

  virtual unsigned int  nFitParams();
  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(paramGamKsKlKK& theParamVal);
  virtual void dumpCurrentResult(std::ostream& os, paramGamKsKlKK& theParamVal, std::string& suffix);

protected:

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);  

  complex<double> etacGammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> eta2225GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> f22340GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> f22300GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> f22010GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> etaToPhiPhiTo4KAmp(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);

private:
  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;

};

#endif
