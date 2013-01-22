// AbsDecayList class definition file. -*- C++ -*-
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
class AbsDecay;

class AbsDecayList {

public:

  AbsDecayList();
  ~AbsDecayList();

  void addDecay(boost::shared_ptr<AbsDecay> theIsoDec);
  boost::shared_ptr<AbsDecay> decay(Particle* mother);
  boost::shared_ptr<AbsDecay> decay(const std::string& name);
  void replaceSuffix(const std::string& oldPart, const std::string& newPart);
  void replaceMassKey(const std::string& oldPart, const std::string& newPart);
  std::vector<boost::shared_ptr<AbsDecay> >& getList() {return _absDecList;}    
protected:  

  std::vector<boost::shared_ptr<AbsDecay> > _absDecList;  
};
