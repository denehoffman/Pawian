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

  std::vector<PbarP*>& allStates(){return _allStates;};
  std::vector<PbarP*>& singletStates(){return _singletStates;};
  std::vector<PbarP*>& tripletM0States(){return _tripletM0States;};
  std::vector<PbarP*>& tripletMp1States(){return _tripletMp1States;};
  std::vector<PbarP*>& tripletMm1States(){return _tripletMm1States;};


 protected:

 private:
  int _jmax;
  std::vector<PbarP*> _allStates;
  std::vector<PbarP*> _singletStates;
  std::vector<PbarP*> _tripletM0States;
  std::vector<PbarP*> _tripletMp1States;
  std::vector<PbarP*> _tripletMm1States;
  jpcRes _pbarJPC;
  jpcRes _pJPC;  
  bool calcJPCs();
};




#endif /* _pbarpStates_H */
