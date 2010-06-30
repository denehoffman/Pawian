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
  pbarpToOmegaPi0States( boost::shared_ptr<pbarpStates> pbarpStates);

  /** Destructor */
  virtual ~pbarpToOmegaPi0States();

  virtual void print(std::ostream& os) const;


 protected:
  virtual bool calcJPCs();

 private:
  jpcRes _omegaJPC;
  jpcRes _piJPC;

  boost::shared_ptr<pbarpStates> _pbarpStates;
};




#endif /* _pbarpToOmegaPi0States_H */
