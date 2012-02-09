#ifndef _OmegaTo3PiLhProdLS_H
#define _OmegaTo3PiLhProdLS_H

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

class OmegaTo3PiLhProdLS : public AbsOmegaPiLhLS{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaTo3PiLhProdLS(boost::shared_ptr<const AbsOmegaPiEventListLS>, boost::shared_ptr<const pbarpToOmegaPi0StatesLS>);
  OmegaTo3PiLhProdLS(boost::shared_ptr<OmegaTo3PiLhProdLS>);

  /** Destructor */
  virtual ~OmegaTo3PiLhProdLS();

  virtual AbsOmegaPiLhLS* clone_() const{
    return new OmegaTo3PiLhProdLS(_omegaPiEventListPtr, _omegaPi0StatesPtr);
        }


  // Getters:
  virtual double calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal);  
  virtual double calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);


  virtual void print(std::ostream& os) const;

protected:
  virtual complex<double> calcCoherentAmp(Spin lamOm, Spin Minit, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData);



private:


};

#endif
