#ifndef _OmegaPiLh_H
#define _OmegaPiLh_H

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

class OmegaPiLh {

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiLh(boost::shared_ptr<const OmegaPiEventList>, boost::shared_ptr<const pbarpToOmegaPi0States>);
  OmegaPiLh(boost::shared_ptr<OmegaPiLh>);

  /** Destructor */
  virtual ~OmegaPiLh();

  OmegaPiLh* clone_() const {
                return new OmegaPiLh(*this);
        }


  // Getters:
  
  double calcLogLh(const OmegaPiData::fitParamVal& theParamVal);
  double calcEvtIntensity(OmegaPiData::OmPiEvtData& theData, const OmegaPiData::fitParamVal& theParamVal);

  boost::shared_ptr<const OmegaPiEventList> getEventList() const {return _omegaPiEventListPtr;}
  boost::shared_ptr<const pbarpToOmegaPi0States> omegaPi0States() const {return _omegaPi0StatesPtr;}
  void print(std::ostream& os) const;

protected:


private:
  boost::shared_ptr<const OmegaPiEventList> _omegaPiEventListPtr;
  boost::shared_ptr<const pbarpToOmegaPi0States> _omegaPi0StatesPtr;

  complex<double> calcCoherentAmp(Spin lamgamma, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData& theData);

};

#endif
