#ifndef _AbsOmegaPiLhLS_H
#define _AbsOmegaPiLhLS_H

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
#include "Examples/pbarpToOmegaPiLS/OmegaPiDataLS.hh"
#include "Examples/pbarpToOmegaPiLS/MinuitstartparamLS.hh"
#include "PwaUtils/DataUtils.hh"

#include <geneva/GConstrainedDoubleObjectCollection.hpp>
#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class AbsOmegaPiEventListLS;
class pbarpToOmegaPi0StatesLS;


class AbsOmegaPiLhLS {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsOmegaPiLhLS(boost::shared_ptr<const AbsOmegaPiEventListLS> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0StatesLS> theStates);
  AbsOmegaPiLhLS(boost::shared_ptr<AbsOmegaPiLhLS> theOmegaPiLhPtr);

  /** Destructor */
  virtual ~AbsOmegaPiLhLS();

  virtual AbsOmegaPiLhLS* clone_() const=0;


  // Getters:
  
  virtual double calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal);
  virtual double calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal)=0;
  boost::shared_ptr<const AbsOmegaPiEventListLS> getEventList() const {return _omegaPiEventListPtr;}
  boost::shared_ptr<const pbarpToOmegaPi0StatesLS> omegaPi0States() const {return _omegaPi0StatesPtr;}

  virtual void getFitParamVal(OmegaPiDataLS::fitParamVal& theParamVal, const std::vector<double>& par) const;
  virtual void setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr );
  virtual void setMnUsrParams(MnUserParameters& upar);
  virtual void setMnUsrParams(MnUserParameters& upar, OmegaPiDataLS::fitParamVal &finalFitParm);
  virtual void setMnUsrParams(MnUserParameters& upar, MinuitStartParamLS &theStartParam);
  virtual void fixMnUsrParams(MnUserParameters& upar, std::vector<std::string>& fixedParams);

  virtual complex<double> spinDensity(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);
  virtual complex<double> spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);


  virtual void dumpCurrentResult(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal) const;
  virtual void printFitParams(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal);
  virtual void print(std::ostream& os) const;

protected:

  unsigned int _globalItCounter;
  boost::shared_ptr<const AbsOmegaPiEventListLS> _omegaPiEventListPtr;
  boost::shared_ptr<const pbarpToOmegaPi0StatesLS> _omegaPi0StatesPtr;
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*> _evtDataVec;
  std::vector<OmegaPiDataLS::OmPiEvtDataLS*> _evtMCVec;

  std::vector< boost::shared_ptr<const JPCLSls> > _allJPCLSlsStates;
  std::vector< boost::shared_ptr<const JPCLSls> > _singlet_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLSls> > _triplet0_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLSls> > _tripletp1_JPCLS_States;
  std::vector< boost::shared_ptr<const JPCLSls> > _tripletm1_JPCLS_States;

  virtual complex<double>  calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLSls> theJPCLSls, pair<double, double> fitVal, OmegaPiDataLS::OmPiEvtDataLS* theData);

  virtual complex<double> calcOmegaProdAmp(Spin Minit, Spin lamomega, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData);
private:

};

#endif
