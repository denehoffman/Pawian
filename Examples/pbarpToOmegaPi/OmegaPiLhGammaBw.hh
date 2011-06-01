#ifndef _OmegaPiLhGammaBw_H
#define _OmegaPiLhGammaBw_H

#include "Examples/pbarpToOmegaPi/OmegaPiLhGamma.hh"



class OmegaPiEventList;
class pbarpToOmegaPi0States;

class OmegaPiLhGammaBw : public OmegaPiLhGamma{

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiLhGammaBw(boost::shared_ptr<const OmegaPiEventList>, boost::shared_ptr<const pbarpToOmegaPi0States>);
  OmegaPiLhGammaBw(boost::shared_ptr<OmegaPiLhGammaBw>);

  /** Destructor */
  virtual ~OmegaPiLhGammaBw();

  virtual AbsOmegaPiLh* clone_() const{
    return new OmegaPiLhGammaBw(_omegaPiEventListPtr, _omegaPi0StatesPtr);
        }


  // Getters:
  virtual double calcLogLh(const OmegaPiData::fitParamVal& theParamVal);  
  virtual double calcEvtIntensity(OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);

  virtual complex<double> spinDensity(Spin M, Spin M_, OmegaPiData::OmPiEvtData* theData, const OmegaPiData::fitParamVal& theParamVal);

  virtual void getFitParamVal(OmegaPiData::fitParamVal& theParamVal, const std::vector<double>& par) const;
  virtual void setGenevaFitParamVal( boost::shared_ptr<Gem::Geneva::GConstrainedDoubleObjectCollection> theGbdc_ptr );
  virtual void setMnUsrParams(MnUserParameters& upar);
  virtual void setMnUsrParams(MnUserParameters& upar, OmegaPiData::fitParamVal &finalFitParm);
  virtual void setMnUsrParams(MnUserParameters& upar, minuitStartParam &theStartParam);

  virtual void dumpCurrentResult(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal) const;
  virtual void printFitParams(std::ostream& os, const OmegaPiData::fitParamVal& fitParmVal);
  virtual void print(std::ostream& os) const;

protected:


private:
  virtual complex<double> calcCoherentAmp(Spin lam, Spin Minit, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm, double fitParmOmMass, double fitParmOmWidth, OmegaPiData::OmPiEvtData* theData);

  complex<double> calcOmegaProdAmp(Spin Minit, Spin lamomega, std::map< boost::shared_ptr<const JPCLS>, pair<double, double>, pawian::Collection::SharedPtrLess >& fitParm,  double fitParmOmMass, double fitParmOmWidth, OmegaPiData::OmPiEvtData* theData);

  complex<double>  calcOmegaProdPartAmp(Spin Minit, Spin lamomega, boost::shared_ptr<const JPCLS> theJPCLS, pair<double, double> fitVal,  double fitParmOmMass, double fitParmOmWidth, OmegaPiData::OmPiEvtData* theData, unsigned int piGamCombi);
};

#endif
