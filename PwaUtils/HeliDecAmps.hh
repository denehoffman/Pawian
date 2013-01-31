// HeliDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/AbsXdecAmp.hh"

class IsobarHeliDecay;
class AbsDecay;

class HeliDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  HeliDecAmps(boost::shared_ptr<IsobarHeliDecay> theDec);
  HeliDecAmps(boost::shared_ptr<AbsDecay> theDec);
  /** Destructor */
  virtual ~HeliDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, Spin lamFs=0);
  virtual complex<double> XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr,
				      EvtData* theData, Spin lamFs);

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual void print(std::ostream& os) const;
  virtual bool checkRecalculation(fitParams& theParamVal);
  // boost::shared_ptr<const jpcRes>& jpcPtr() {return _JPCPtr;}
  std::vector< boost::shared_ptr<const JPClamlam> >& jpclamlamVec() {return _JPClamlams;}
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  std::vector< boost::shared_ptr<const JPClamlam> > _JPClamlams;
  double _factorMag;
  double _parityFactor;
  std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > _currentParamMagLamLams;
  std::map< boost::shared_ptr<const JPClamlam>, double, pawian::Collection::SharedPtrLess > _currentParamPhiLamLams;
private:




};



