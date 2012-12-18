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

class WignerDUtil {

public:

  WignerDUtil(std::vector<Particle*>& mother, std::vector<Particle*>& daughter);
  ~WignerDUtil();

  std::string name() {return _name;}
    
protected:  
  std::string _name;
private:

};
