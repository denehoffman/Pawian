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
  JpsiGamKsKlKKProdLh(boost::shared_ptr<const JpsiGamKsKlKKEventList>, const std::map<const std::string, bool>& hypMap);
  JpsiGamKsKlKKProdLh(boost::shared_ptr<AbsJpsiGamKsKlKKLh>, const std::map<const std::string, bool>& hypMap);

  /** Destructor */
  virtual ~JpsiGamKsKlKKProdLh();

  virtual AbsJpsiGamKsKlKKLh* clone_() const {
    return new JpsiGamKsKlKKProdLh(_JpsiGamKsKlKKEvtListPtr, _hypMap);
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

  complex<double> etaGammaAmp(Spin Minit, Spin Metac, Spin Mgamma, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampParam, pair<double, double>& massParam);

  complex<double> f0GammaAmp(Spin Minit, Spin Mgamma, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampf0Prod, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampf0Dec, pair<double, double>& massParam);



  complex<double> f22340GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> f22300GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> f22010GammaCoherentAmp(Spin Minit, Spin Metac, Spin Mgamma, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> etaToPhiPhiTo4KAmp(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData);
  complex<double> f0ToPhiPhiTo4KAmp(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess > ampf0Dec);


  bool _eta2225Hyp;
  bool _f02020Hyp;
  std::map<const std::string, bool> _hypMap;
private:
  unsigned int _nFitParams;
  std::vector<unsigned int> _ampVec;
  std::vector<unsigned int> _massVec;

};

#endif
