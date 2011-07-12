#ifndef _OmegaPiLhPi0GammaLS_H
#define _OmegaPiLhPi0GammaLS_H

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

class OmegaPiLhPi0GammaLS : public AbsOmegaPiLhLS{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiLhPi0GammaLS(boost::shared_ptr<const AbsOmegaPiEventListLS>, boost::shared_ptr<const pbarpToOmegaPi0StatesLS>);
  OmegaPiLhPi0GammaLS(boost::shared_ptr<OmegaPiLhPi0GammaLS>);

  /** Destructor */
  virtual ~OmegaPiLhPi0GammaLS();

  virtual AbsOmegaPiLhLS* clone_() const{
    return new OmegaPiLhPi0GammaLS(_omegaPiEventListPtr, _omegaPi0StatesPtr);
        }


  // Getters:
  virtual double calcLogLh(const OmegaPiDataLS::fitParamVal& theParamVal);  
  virtual double calcEvtIntensity(OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);

  virtual complex<double> spinDensity(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);
  virtual complex<double> spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiDataLS::OmPiEvtDataLS* theData, const OmegaPiDataLS::fitParamVal& theParamVal);

  virtual void getFitParamVal(OmegaPiDataLS::fitParamVal& theParamVal, const std::vector<double>& par) const;
  virtual void setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr );
  virtual void setMnUsrParams(MnUserParameters& upar);
  virtual void setMnUsrParams(MnUserParameters& upar, OmegaPiDataLS::fitParamVal &fitParam);
  virtual void setMnUsrParams(MnUserParameters& upar, MinuitStartParamLS &theStartParam);

  virtual void dumpCurrentResult(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal) const;
  virtual void printFitParams(std::ostream& os, const OmegaPiDataLS::fitParamVal& fitParmVal);
  virtual void print(std::ostream& os) const;

protected:
  virtual complex<double> calcCoherentAmp(Spin lamgamma, Spin Minit, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData);

protected:


private:
  complex<double>  calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLSls> theJPCLSls, pair<double, double> fitVal, OmegaPiDataLS::OmPiEvtDataLS* theData);

  complex<double> calcOmegaProdAmp(Spin Minit, Spin lamomega, const OmegaPiDataLS::fitParamVal& theParamVal, std::vector< boost::shared_ptr<const JPCLSls> >& theJPCLSlsStates, OmegaPiDataLS::OmPiEvtDataLS* theData);


};

#endif
