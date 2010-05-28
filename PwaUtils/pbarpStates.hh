#ifndef _pbarpStates_H
#define _pbarpStates_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class pbarpStates {

public:

  /// Default Constructor
  pbarpStates();
  pbarpStates(int jmax);

  /** Destructor */
  virtual ~pbarpStates();

  void print(std::ostream& os) const;

//   std::vector<jpcRes*> allStates();
//   std::vector<jpcRes*> lrange(int lmin; int lmax);
//   std::vector<jpcRes*> jrange(Spin jmin; Spin jmax);

 protected:

 private:
  int _jmax;
//   std::vector<jpcRes*> _theStates;
  std::vector<PbarP*> _theStates;
  jpcRes _pbarJPC;
  jpcRes _pJPC;  
  bool calcJPCs();
};




#endif /* _pbarpStates_H */
