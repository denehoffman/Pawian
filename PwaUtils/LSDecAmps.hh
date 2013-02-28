// LSDecAmps class definition file. -*- C++ -*-
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

class IsobarLSDecay;
class AbsDecay;

class LSDecAmps : public AbsXdecAmp{

public:

  // create/copy/destroy:

  ///Constructor 
  LSDecAmps(boost::shared_ptr<IsobarLSDecay> theDec);
  LSDecAmps(boost::shared_ptr<AbsDecay> theDec);
  /** Destructor */
  virtual ~LSDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, Spin lamFs);
  virtual complex<double> XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr,
				      EvtData* theData, Spin lamFs);

  virtual void print(std::ostream& os) const;
  std::vector< boost::shared_ptr<const JPCLS> >& jpclsVec() {return _JPCLSs;}

  virtual void getDefaultParams(fitParams& fitVal, fitParams& fitErr);
  virtual bool checkRecalculation(fitParams& theParamVal);
  virtual void updateFitParams(fitParams& theParamVal);

protected:
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSs;
  double _factorMag;
  double _parityFactor;
  std::map< boost::shared_ptr<const JPCLS>, map<Spin,map<Spin, double > > > _cgPreFactor;

  void  fillCgPreFactor();
  virtual complex<double> lsLoop(Spin lamX, EvtData* theData, Spin lam1Min, Spin lam1Max, Spin lam2Min, Spin lam2Max, bool withDecs, Spin lamFs=0 ); 
private:




};



