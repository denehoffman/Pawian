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
#include "PwaUtils/DataUtils.hh"

// using namespace std;


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
  
  double calcLogLh(const OmegaPiData::fitParamVal& theParamVal);
  virtual double calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)=0;
  virtual complex<double> spinDensity(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)=0;
  virtual complex<double> spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal)=0;
  boost::shared_ptr<const OmegaPiEventList> getEventList() const {return _omegaPiEventListPtr;}
  boost::shared_ptr<const pbarpToOmegaPi0States> omegaPi0States() const {return _omegaPi0StatesPtr;}
  virtual void print(std::ostream& os) const;

protected:

  boost::shared_ptr<const OmegaPiEventList> _omegaPiEventListPtr;
  boost::shared_ptr<const pbarpToOmegaPi0States> _omegaPi0StatesPtr;
  std::vector<OmegaPiData::OmPiEvtData*> _evtDataVec;
  std::vector<OmegaPiData::OmPiEvtData*> _evtMCVec;

  virtual complex<double> calcCoherentAmp(Spin lam, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData) = 0;

private:

};

#endif
