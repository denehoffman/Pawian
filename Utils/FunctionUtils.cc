#include "Utils/FunctionUtils.hh"
#include "Particle/Particle.hh"


namespace FunctionUtils{


  std::string particleListName(std::vector<Particle*>& particleVec){
  std::string result;

  std::vector<Particle*>::iterator it;
  for(it=particleVec.begin(); it!=particleVec.end(); ++it){
    result+=(*it)->name();
  }

  return result;
  }

}
