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
  static IsobarDecayList* instance();

  ~IsobarDecayList();

  void addDecay(boost::shared_ptr<IsobarDecay> theIsoDec);
  boost::shared_ptr<IsobarDecay> decay(Particle* mother);
 
protected:  

  IsobarDecayList(){;}
  static IsobarDecayList* _instance;
  std::map<std::string, boost::shared_ptr<IsobarDecay> > _isoDecList;  
};
