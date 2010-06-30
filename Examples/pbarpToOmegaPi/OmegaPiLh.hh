#ifndef _OmegaPiLh_H
#define _OmegaPiLh_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>
#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "Examples/pbarpToOmegaPi/OmegaPiData.hh"
#include "PwaUtils/DataUtils.hh"

// using namespace std;


class OmegaPiEventList;


class OmegaPiLh {

public:

  // create/copy/destroy:

  ///Constructor 
  OmegaPiLh(boost::shared_ptr<const OmegaPiEventList>, int);
  OmegaPiLh(boost::shared_ptr<OmegaPiLh>);

  /** Destructor */
  virtual ~OmegaPiLh();

  OmegaPiLh* clone_() const {
                return new OmegaPiLh(*this);
        }


  // Getters:
  
  //   double calcLogLh(const OmegaPiData::fitParamVal& theParamVal);
  //   double calcEvtIntensity(const OmegaPiData::OmPiEvtData& theEvtData, const fitParamVal& theParamVal);
  //   bool setFitParamVal(OmegaPiData::fitParamVal& fitParamVal, const std::vector<double>& par);
  boost::shared_ptr<const OmegaPiEventList> getEventList() const {return _omegaPiEventListPtr;}
  std::vector< boost::shared_ptr<const JPCLS> > omegaProdStates() const {return _JPCLSomegaProd;}
  int jmax() const {return _jmax;}
  void print(std::ostream& os) const;

protected:


private:
  boost::shared_ptr<const OmegaPiEventList> _omegaPiEventListPtr;
  int _jmax;

  std::vector< boost::shared_ptr<const jpcRes> > _pbarpJPC;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProd;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpSinglett;

  void init();

};

#endif
