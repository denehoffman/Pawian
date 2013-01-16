// IsobarDecayList class definition file. -*- C++ -*-
// Copyright 2012 Bertram Kopf

#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <boost/shared_ptr.hpp>

#include "PwaUtils/DataUtils.hh"

class Particle;
class IsobarDecay;

class IsobarDecayList {

public:

  IsobarDecayList();
  ~IsobarDecayList();

  void addDecay(boost::shared_ptr<IsobarDecay> theIsoDec);
  boost::shared_ptr<IsobarDecay> decay(Particle* mother);
  boost::shared_ptr<IsobarDecay> decay(const std::string& name);
  void replaceSuffix(const std::string& oldPart, const std::string& newPart);
  void replaceMassKey(const std::string& oldPart, const std::string& newPart);
  std::vector<boost::shared_ptr<IsobarDecay> >& getList() {return _isoDecList;}    
protected:  

  std::vector<boost::shared_ptr<IsobarDecay> > _isoDecList;  
};
