#ifndef _AbsOmegaPiLh_H
#define _AbsOmegaPiLh_H

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
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "Examples/pbarpToOmegaPi/minuitstartparam.hh"
#include "PwaUtils/DataUtils.hh"

#include <geneva/GConstrainedDoubleObjectCollection.hpp>
#include "Minuit2/MnUserParameters.h"

// using namespace std;
using namespace ROOT::Minuit2;

class OmegaPiEventList;
class pbarpToOmegaPi0States;


class AbsOmegaPiLh {

public:

  // create/copy/destroy:

  ///Constructor 
  AbsOmegaPiLh(boost::shared_ptr<const OmegaPiEventList> theEvtList, boost::shared_ptr<const pbarpToOmegaPi0States> theStates);
  AbsOmegaPiLh(boost::shared_ptr<AbsOmegaPiLh> theOmegaPiLhPtr);

  /** Destructor */
  virtual ~AbsOmegaPiLh();

  virtual AbsOmegaPiLh* clone_() const=0;


  // Getters:
  
  virtual double calcLogLh(const OmegaPiData::fitParamVal& theParamVal);
  virtual double calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)=0;
  virtual complex<double> spinDensity(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)=0;
  virtual complex<double> spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)=0;
  boost::shared_ptr<const OmegaPiEventList> getEventList() const {return _omegaPiEventListPtr;}
  boost::shared_ptr<const pbarpToOmegaPi0States> omegaPi0States() const {return _omegaPi0StatesPtr;}

  virtual void getFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const=0;
  virtual void setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr )=0;
  virtual void setMnUsrParams(MnUserParameters& upar)=0;
  virtual void setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm)=0;
  virtual void setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam)=0;

  virtual void dumpCurrentResult(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal) const =0;
  virtual void printFitParams(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal)=0;
  virtual void print(std::ostream& os) const;

protected:
  unsigned int _globalItCounter;
  boost::shared_ptr<const OmegaPiEventList> _omegaPiEventListPtr;
  boost::shared_ptr<const pbarpToOmegaPi0States> _omegaPi0StatesPtr;
  std::vector<OmegaPiData::OmPiEvtData*> _evtDataVec;
  std::vector<OmegaPiData::OmPiEvtData*> _evtMCVec;


private:

};

#endif
