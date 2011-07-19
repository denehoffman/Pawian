#ifndef _JPCdecays_H
#define _JPCdecays_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <cassert>

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

class jpcRes;

class JPCdecays {

public:

  // create/copy/destroy:

  ///Constructor 
  JPCdecays(boost::shared_ptr<const jpcRes> motherJPC, boost::shared_ptr<const jpcRes> daughter1JPC, boost::shared_ptr<const jpcRes> daughter2JPC, bool sameDaughters=false);


  /** Destructor */
  virtual ~JPCdecays();

  // Getters:
  void print(std::ostream& os) const; 

protected:


private:
  boost::shared_ptr<const jpcRes> _motherJPC;
  boost::shared_ptr<const jpcRes> _daughter1JPC;
  boost::shared_ptr<const jpcRes> _daughter2JPC;
  std::vector<LS> _allLSs;
  std::vector<LS> _allowedCParityLSs;
  std::vector<LS> _forbiddenCParityLSs;
  bool _sameDaughters;
};

#endif
