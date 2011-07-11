#ifndef _pbarpToOmegaPi0States_H
#define _pbarpToOmegaPi0States_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/AbsStates.hh"



class pbarpStates;

class pbarpToOmegaPi0States : public AbsStates {

public:

  /// Default Constructor  
  pbarpToOmegaPi0States(boost::shared_ptr<pbarpStates> pbarpStates);

  /** Destructor */
  virtual ~pbarpToOmegaPi0States();

  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpSingletStates() const {return _pbarpSingletToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpTripletM0States() const {return _pbarpTripletM0ToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const jpcRes> > pbarpTripletM1States() const {return _pbarpTripletM1ToOmegaPi;};
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsSinglet() const {return _JPCLSomegaProdSinglet;};
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet0() const {return _JPCLSomegaProdTripletM0;};
  virtual std::vector< boost::shared_ptr<const JPCLS> > jpclsTriplet1() const {return _JPCLSomegaProdTripletM1;};

  virtual void print(std::ostream& os) const;


 protected:
  virtual bool calcStates();
  std::vector< boost::shared_ptr<const JPCLSM> > _allStates;
  std::vector< boost::shared_ptr<const JPCSM> > _allJPCSM;
  std::vector< boost::shared_ptr<const JPCLS> > _allJPCLS;
  std::vector< boost::shared_ptr<const jpcRes> >  _jpcStates;

 private:
  jpcRes _omegaJPC;
  jpcRes _piJPC;

  boost::shared_ptr<pbarpStates> _pbarpStatesAll;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpSingletToOmegaPi;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM0ToOmegaPi;
  std::vector< boost::shared_ptr<const jpcRes> > _pbarpTripletM1ToOmegaPi;
  std::vector< boost::shared_ptr<const JPCLS> >  _JPCLSomegaProdSinglet;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProdTripletM0;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLSomegaProdTripletM1;
  std::vector< boost::shared_ptr<const JPCLS> > extractJPCLSStates(std::vector< boost::shared_ptr<const jpcRes> > , std::vector< boost::shared_ptr<const JPCLS> > ) const;
};




#endif /* _pbarpToOmegaPi0States_H */
