#ifndef _OmegaTo3PiLhPi0GammaLS_H
#define _OmegaTo3PiLhPi0GammaLS_H

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
#include "Examples/pbarpToOmegaPiLS/AbsOmegaPiLhLS.hh"


#include "PwaUtils/DataUtils.hh"

// using namespace std;


class AbsOmegaPiEventListLS;
class pbarpToOmegaPi0StatesLS;

class OmegaTo3PiLhPi0GammaLS : public AbsOmegaPiLhLS{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaTo3PiLhPi0GammaLS(boost::shared_ptr<const AbsOmegaPiEventListLS>, boost::shared_ptr<const pbarpToOmegaPi0StatesLS>);
  OmegaTo3PiLhPi0GammaLS(boost::shared_ptr<OmegaTo3PiLhPi0GammaLS>);

  /** Destructor */
  virtual ~OmegaTo3PiLhPi0GammaLS();

  virtual AbsOmegaPiLhLS* clone_() const{
    return new OmegaTo3PiLhPi0GammaLS(_omegaPiEventListPtr, _omegaPi0StatesPtr);
        }


  // Getters:
  virtual double calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal);  
  virtual double calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);


  virtual void print(std::ostream& os) const;

protected:
  virtual complex<double> calcCoherentAmp(Spin, Spin, const OmegaPiDataLS::fitParamVal& , std::vector< boost::shared_ptr<const JPCLSls> >& , OmegaPiDataLS::OmPiEvtDataLS*);



private:


};

#endif
