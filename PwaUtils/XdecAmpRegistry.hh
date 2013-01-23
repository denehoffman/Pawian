// DecAmpRegistry class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <string>

#include <cassert>
#include <boost/shared_ptr.hpp>

class AbsXdecAmp;
class AbsDecay;

class XdecAmpRegistry{

public:

  /** Destructor */
  virtual ~XdecAmpRegistry();

  static XdecAmpRegistry* instance();
  boost::shared_ptr<AbsXdecAmp> getXdecAmp(boost::shared_ptr<AbsDecay>);

protected:
 ///Constructor 
  XdecAmpRegistry();
  static XdecAmpRegistry* _instance;

private:

  std::map<std::string, boost::shared_ptr<AbsXdecAmp> > _xDecAmpMap;



};



