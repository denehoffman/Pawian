#ifndef _etacToapi0lh_H
#define _etacToapi0lh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Data.hh"

// using namespace std;


class EtacToapi0EventList;


class EtacToapi0Lh {

public:

  // create/copy/destroy:

  ///Constructor 
  EtacToapi0Lh(boost::shared_ptr<const EtacToapi0EventList>);
  EtacToapi0Lh(boost::shared_ptr<EtacToapi0Lh>);

  /** Destructor */
  virtual ~EtacToapi0Lh();

  EtacToapi0Lh* clone_() const {
                return new EtacToapi0Lh(*this);
        }
  // Getters:

  double calcLogLh(const fitParamVal& theParamVal);
  double calcEvtIntensity(const evt4Vec& theEvtVec, const fitParamVal& theParamVal);
  bool setFitParamVal(fitParamVal& fitParamVal, const std::vector<double>& par);
  boost::shared_ptr<const EtacToapi0EventList> getEventList() const {return _etacToapi0EventListPtr;}

protected:


private:
  boost::shared_ptr<const EtacToapi0EventList> _etacToapi0EventListPtr;

};

#endif
