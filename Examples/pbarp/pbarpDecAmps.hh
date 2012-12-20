#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"

class IsobarDecay;

class pbarpDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  pbarpDecAmps(boost::shared_ptr<IsobarDecay> theDec);

  /** Destructor */
  virtual ~pbarpDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtDataNew* theData);
  virtual complex<double> daughterAmp(Spin lam1, Spin lam2, EvtDataNew* theData);                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParamsNew& theParamVal);
  boost::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsVec() {return _JPCLSs;}
  void updateFitParams(fitParamsNew& theParamVal);

protected:
  boost::shared_ptr<IsobarDecay> _decay; 
  boost::shared_ptr<const jpcRes> _JPCPtr;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSs;
  boost::shared_ptr<pbarpDecAmps> _decAmpDaughter1;
  boost::shared_ptr<pbarpDecAmps> _decAmpDaughter2;
  Spin _Jdaughter1;
  Spin _Jdaughter2;
  std::string _key;
  std::string _massKey;
  const std::string _wignerDKey;
  bool _daughter1IsStable;
  bool _daughter2IsStable;
  bool _withDyn;
  virtual void initialize();

private:




};



