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
                                         
  virtual void getDefaultParams(fitParamsNew& fitVal, fitParamsNew& fitErr);
  virtual void print(std::ostream& os) const;
  virtual void checkRecalculation(fitParamsNew& theParamVal);

  void updateFitParams(fitParamsNew& theParamVal);
protected:
  boost::shared_ptr<IsobarDecay> _decay;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSs;
  boost::shared_ptr<pbarpDecAmps> _decAmpDaughter1;
  boost::shared_ptr<pbarpDecAmps> _decAmpDaughter2;
  std::string _key;
  bool _daughter1IsStable;
  bool _daughter2IsStable;
  virtual void initialize();

private:




};



