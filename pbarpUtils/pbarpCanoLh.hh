// pbarpCanoLh class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "pbarpUtils/pbarpBaseLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class pbarpReaction;
class LSDecAmps;

class pbarpCanoLh : public pbarpBaseLh {

public:
  pbarpCanoLh(boost::shared_ptr<const EvtDataBaseList>);

  virtual ~pbarpCanoLh();

  virtual AbsLh* clone_() const {
    return new  pbarpCanoLh(_evtListPtr);
  }
  virtual void print(std::ostream& os) const;


protected:


private:

  void initialize();
};
