#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include <boost/shared_ptr.hpp>

// #include <TSystem.h>
#include "qft++/topincludes/relativistic-quantum-mechanics.hh"
#include "PwaUtils/DataUtils.hh"


class D0ToKsPipPimStates {

public:

  ///Constructor 
  D0ToKsPipPimStates();

  /** Destructor */
  virtual ~D0ToKsPipPimStates();

  std::vector< boost::shared_ptr<const JPCLS> > D0Tof0KStates() const {return _JPCLS_D0ToF0K;}
  std::vector< boost::shared_ptr<const JPCLS> > D0ToK0KStates() const {return _JPCLS_D0ToK0K;}
 
  virtual void print(std::ostream& os) const; 

protected:
  virtual void fillJPCLS(boost::shared_ptr<jpcRes>, boost::shared_ptr<jpcRes>, 
                         boost::shared_ptr<jpcRes>, 
                         std::vector< boost::shared_ptr<const JPCLS> >& );
virtual void printDecayJPCLS(std::ostream& os,
                               std::vector< boost::shared_ptr<const JPCLS > > theJPCLS) const;
private:
  boost::shared_ptr<jpcRes> _D0JPC;
  boost::shared_ptr<jpcRes> _f0JPC;
  boost::shared_ptr<jpcRes> _piJPC;
  boost::shared_ptr<jpcRes> _kJPC;
  boost::shared_ptr<jpcRes> _k0JPC;

  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_D0ToF0K;
  std::vector< boost::shared_ptr<const JPCLS> > _JPCLS_D0ToK0K;


};

