#ifndef _pbarpStates_H
#define _pbarpStates_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/AbsStates.hh"
#include "PwaUtils/DataUtils.hh"


class pbarpStates : public AbsStates {

public:

  /// Default Constructor
  pbarpStates();
  pbarpStates(int jmax);

  /** Destructor */
  virtual ~pbarpStates();

  virtual void print(std::ostream& os) const;

  virtual std::vector< boost::shared_ptr<const JPCLSM> > singletStates() {return _singletStates;};
  virtual std::vector< boost::shared_ptr<const JPCLSM> > tripletM0States() {return _tripletM0States;};
  virtual std::vector< boost::shared_ptr<const JPCLSM> > tripletMp1States() {return _tripletMp1States;};
  virtual std::vector< boost::shared_ptr<const JPCLSM> > tripletMm1States() {return _tripletMm1States;};

 protected:
  virtual bool calcJPCs();

 private:
  int _jmax;
  jpcRes _pbarJPC;
  jpcRes _pJPC;

  std::vector< boost::shared_ptr<const JPCLSM> > _singletStates;
  std::vector< boost::shared_ptr<const JPCLSM> > _tripletM0States;
  std::vector< boost::shared_ptr<const JPCLSM> > _tripletMp1States;
  std::vector< boost::shared_ptr<const JPCLSM> > _tripletMm1States;  
};




#endif /* _pbarpStates_H */
