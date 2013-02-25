// DynRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

class AbsDecay;
class AbsDynamics;
class Particle;

class DynRegistry{

public:

  /** Destructor */
  virtual ~DynRegistry();

  static DynRegistry* instance();
  boost::shared_ptr<AbsDynamics> getDynamics(boost::shared_ptr<AbsDecay> theDec);

protected:
 ///Constructor 
  DynRegistry();
  static DynRegistry* _instance;

private:

  std::map<std::string, boost::shared_ptr<AbsDynamics> > _dynMap;


};



