#ifndef _pbarpToOmegaPi0States_H
#define _pbarpToOmegaPi0States_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"

class pbarpStates;

class pbarpToOmegaPi0States : public AbsStates {

public:

  /// Default Constructor
  pbarpToOmegaPi0States(boost::shared_ptr<pbarpStates> pbarpStates);

  /** Destructor */
  virtual ~pbarpToOmegaPi0States();

  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpSingletStates() const {return _pbarpSingletToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpTripletM0States() const {return _pbarpTripletM0ToOmegaPi;};
  virtual void print(std::ostream& os) const;


 protected:
  virtual bool calcJPCs();

 private:
  jpcRes _omegaJPC;
  jpcRes _piJPC;

  boost::shared_ptr<pbarpStates> _pbarpStatesAll;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpSingletToOmegaPi;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM0ToOmegaPi;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM1ToOmegaPi;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProd;
};




#endif /* _pbarpToOmegaPi0States_H */
