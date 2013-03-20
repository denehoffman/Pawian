// epemBaseLh class definition file. -*- C++ -*-
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

#include "PwaUtils/AbsLh.hh"
#include "PwaUtils/DataUtils.hh"
#include "Minuit2/MnUserParameters.h"

class AbsXdecAmp;
class epemReaction;
class LSDecAmps;

class epemBaseLh : public AbsLh {

public:
  epemBaseLh(boost::shared_ptr<const EvtDataBaseList>);
  
  virtual ~epemBaseLh();
  
  virtual AbsLh* clone_() const{
    return new  epemBaseLh(_evtListPtr);
  }

  virtual double calcEvtIntensity( EvtData* theData, fitParams& theParamVal);
  virtual complex<double> calcProdPartAmp(Spin lamX, Spin lamDec, std::string nameDec, EvtData* theData, 
					  std::map <boost::shared_ptr<const JPCLS>,
					  std::vector< boost::shared_ptr<AbsXdecAmp> >,
					  pawian::Collection::SharedPtrLess > pbarpAmps);

  virtual complex<double> calcSpinDensity(Spin M1, Spin M2, std::string& nameDec, EvtData* theData);

  virtual void print(std::ostream& os) const;

  
protected:
  boost::shared_ptr<epemReaction> _epemReactionPtr;

  int _highestJFsp;
  bool _isHighestJaPhoton;
  
private:

  void initialize();
};
