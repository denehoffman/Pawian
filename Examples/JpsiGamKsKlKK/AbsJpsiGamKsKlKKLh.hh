#ifndef _AbsJpsiGamKsKlKKLh_H
#define _AbsJpsiGamKsKlKKLh_H

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
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKData.hh"
#include "Examples/JpsiGamKsKlKK/JpsiGamKsKlKKFitParams.hh"
#include "PwaUtils/DataUtils.hh"

#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class JpsiGamKsKlKKEventList;

class AbsJpsiGamKsKlKKLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsJpsiGamKsKlKKLh(boost::shared_ptr<const JpsiGamKsKlKKEventList>);
  AbsJpsiGamKsKlKKLh(boost::shared_ptr<AbsJpsiGamKsKlKKLh>);

  /** Destructor */
  virtual ~AbsJpsiGamKsKlKKLh();

  virtual AbsJpsiGamKsKlKKLh* clone_() const =0;


  // Getters:
  
  virtual double calcLogLh(const paramGamKsKlKK& theParamVal);
  virtual double calcEvtIntensity(JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData, const paramGamKsKlKK& theParamVal);

  virtual boost::shared_ptr<const JpsiGamKsKlKKEventList> getEventList() const {return _JpsiGamKsKlKKEvtListPtr;}

  virtual void setMnUsrParams(MnUserParameters& upar, paramGamKsKlKK& startVal,  paramGamKsKlKK& errVal)=0;

  virtual int setFitParamVal(paramGamKsKlKK& theParamVal, const std::vector<double>& par) =0;
  virtual unsigned int nFitParams() =0;

  virtual void print(std::ostream& os) const;
  virtual void printCurrentFitResult(paramGamKsKlKK& theParamVal) =0;
  virtual void dumpCurrentResult(std::ostream& os, paramGamKsKlKK& theParamVal, std::string& suffix)=0;
protected:

  boost::shared_ptr<const JpsiGamKsKlKKEventList> _JpsiGamKsKlKKEvtListPtr;
  JpsiGamKsKlKKFitParams _fitparamsGamKsKlKK; 
 
  std::vector<JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData*> _evtDataVec;
  std::vector<JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData*> _evtMCVec;

  virtual complex<double> calcCoherentAmp(Spin Minit, Spin lamGam, const paramGamKsKlKK& theParamVal, JpsiGamKsKlKKData::JpsiGamKsKlKKEvtData* theData)=0;
private:


};

#endif
