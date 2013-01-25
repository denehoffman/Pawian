// LSOmegaTo3PiDecAmps class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/LSDecAmps.hh"

class OmegaTo3PiLSDecay;

class LSOmegaTo3PiDecAmps : public LSDecAmps{

public:

  // create/copy/destroy:

  ///Constructor 
  LSOmegaTo3PiDecAmps(boost::shared_ptr<OmegaTo3PiLSDecay> theDec);

  /** Destructor */
  virtual ~LSOmegaTo3PiDecAmps();


  // Getters:
  
  virtual complex<double> XdecAmp(Spin lamX, EvtData* theData, Spin lamFs=0);
  virtual complex<double> XdecPartAmp(Spin lamX, Spin lamDec, short fixDaughterNr,
				      EvtData* theData, Spin lamFs);

  virtual void print(std::ostream& os) const;

protected:

  const std::string _lambdaDecKey;

private:




};



