#ifndef _OmegaPiLhOmega_H
#define _OmegaPiLhOmega_H

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
#include "Examples/pbarpToOmegaPi/AbsOmegaPiLh.hh"
#include "PwaUtils/DataUtils.hh"

// using namespace std;


class OmegaPiEventList;
class pbarpToOmegaPi0States;

class OmegaPiLhOmega : public AbsOmegaPiLh{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiLhOmega(boost::shared_ptr<const OmegaPiEventList>, boost::shared_ptr<const pbarpToOmegaPi0States>);
  OmegaPiLhOmega(boost::shared_ptr<OmegaPiLhOmega>);

  /** Destructor */
  virtual ~OmegaPiLhOmega();

  virtual AbsOmegaPiLh* clone_() const{
    return new OmegaPiLhOmega(_omegaPiEventListPtr, _omegaPi0StatesPtr);
        }


  // Getters:
  
  virtual double calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);
  virtual void print(std::ostream& os) const;
protected:
  virtual complex<double> calcCoherentAmp(Spin lamOmega, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData);

private:

};

#endif
