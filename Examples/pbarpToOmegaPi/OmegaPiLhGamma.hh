#ifndef _OmegaPiLhGamma_H
#define _OmegaPiLhGamma_H

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

class OmegaPiLhGamma : public AbsOmegaPiLh{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiLhGamma(boost::shared_ptr<const OmegaPiEventList>, boost::shared_ptr<const pbarpToOmegaPi0States>);
  OmegaPiLhGamma(boost::shared_ptr<OmegaPiLhGamma>);

  /** Destructor */
  virtual ~OmegaPiLhGamma();

  virtual AbsOmegaPiLh* clone_() const{
    return new OmegaPiLhGamma(_omegaPiEventListPtr, _omegaPi0StatesPtr);
        }


  // Getters:
  virtual double calcLogLh(const OmegaPiData::fitParamVal& theParamVal);  
  virtual double calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);

  virtual complex<double> spinDensity(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);
  virtual complex<double> spinDensityOmegaFrame(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);

  virtual void getFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const;
  virtual void setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr );
  virtual void setMnUsrParams(MnUserParameters& upar);
  virtual void setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm);
  virtual void setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam);

  virtual void dumpCurrentResult(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal) const;
  virtual void printFitParams(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal);
  virtual void print(std::ostream& os) const;

protected:
  virtual complex<double> calcCoherentAmp(Spin lam, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData);

protected:


private:
  complex<double>  calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLS> theJPCLS, pair<double, double> fitVal, OmegaPiData::OmPiEvtData* theData);

  complex<double> calcOmegaProdAmp(Spin Minit, Spin lamomega, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, OmegaPiData::OmPiEvtData* theData);

  void FillUpars(MnUserParameters& upar, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParmS, const string &theSuffix);

  void FillUpars(MnUserParameters& upar, minuitStartParam &theStartParam, const std::vector< boost::shared_ptr<const JPCLS> > &theJPCLS, const string &theSuffix);
};

#endif
